classdef CorrelationMap < handle
	
	% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
	% Copyright 2009-2014 Adriano Vilela Barbosa

	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
	% 	VERY IMPORTANT: THE CORRELATION MAP IS STRUCTURED IN SUCH A WAY THAT:
	% 	- A POSITIVE OFFSET MEANS THAT SIGNAL 1 IS LEADING (SIGNAL 2 LAGGING)
	% 	- A NEGATIVE OFFSET MEANS THAT SIGNAL 2 IS LEADING (SIGNAL 1 LAGGING)
	
	% Dependent = true    -> property computed only when asked
	% SetAccess = private -> read-only property
	% GetAccess = private -> property visible only within the class
	
	% Public properties
	properties
		name               = '';
		% input_signal_names = {'' ''};
		start_sample       = 1;
		eta                = 0.1;
		delta              = 0;
		filter_type        = 'ds';
	end
	
	% Private, non-dependent properties
	properties(SetAccess = private)
		signal_1    = [];
		signal_2    = [];
		rate        = [];
		n_samples   = [];
		corr_map    = [];
		corr_matrix = [];
		corr_signal = [];
	end
	
	% Private, dependent properties
	% Important: strictly speaking, the properties 'corr_map', 'corr_matrix'
	% and 'corr_signal' are dependent, in the sense that their values depend
	% on the properties 'eta', 'delta' and 'filter_type'. However, declaring
	% those properties as dependent would mean that they are not stored with
	% the object but, instead, they are computed every time they are queried.
	% Since the computation of the correlation matrix can be time consuming,
	% this becomes a problem. A better approach, which we use here, is to
	% use events and listeners. Basically, we trigger an event every time
	% any of the properties 'eta', 'delta' and 'filter_type' changes. We
	% then create a listener for these events, and attach to it a callback
	% function that updates the correlation matrix accordingly.
	properties(Dependent = true)
		time_offsets = [];
		time_vector  = [];
		time_limits  = [];
	end
	
	% 	Property 'corr_matrix' is the 2D correlation matrix, as returned by
	% 	the function instantaneous_correlation(), with the correlation signals
	% 	along its diagonals. In turn, property 'corr_map' is the correlation
	% 	map obtained by rotating 'corr_matrix' 90 degrees counter-clockwise
	% 	and, therefore, contains the correlation signals along its rows.
	% 	The latter is the property actually used when plotting.
	
	events
		outdated_corr_matrix
	end
	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
	methods
		
		% --------------------------------------------------------------------------------------------- %
		
		function corr_map = CorrelationMap(signal_1,signal_2,filter_type,eta,delta)
			
			% This is the constructor function
			
			% Initialize the input signals
			set_signals(corr_map,signal_1,signal_2);
			
			% Set the property 'eta', if it has been provided
			if (exist('eta','var') & ~isempty(eta))
				corr_map.eta = eta;
			end
			
			% Set the property 'delta', if it has been provided
			if (exist('delta','var') & ~isempty(delta))
				corr_map.delta = delta;
			end
			
			% Set the property 'filter_type', if it has been provided
			if (exist('filter_type','var') & ~isempty(filter_type))
				corr_map.filter_type = filter_type;
			end
			
			% Define event listeners. We could use property events, but class
			% events are much more flexible
			addlistener(corr_map,'outdated_corr_matrix',@update_corr_matrix);
			
			% At this point, the correlation matrix is empty, and we then call
			% the handler for the event 'outdated_corr_matrix' to fill it
			update_corr_matrix(corr_map);
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_vector = get.time_vector(corr_map)
			% The time vector (signal_1 is used as the reference signal)
			time_vector = corr_map.signal_1.time_vector;
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_limits = get.time_limits(corr_map)
			time_limits = corr_map.signal_1.time_limits;
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_offsets = get.time_offsets(corr_map)
			time_offsets = [-corr_map.delta:corr_map.delta]/corr_map.rate;
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function set.eta(corr_map,eta)
			
			if (corr_map.eta ~= eta)
				corr_map.eta = eta;
				notify(corr_map,'outdated_corr_matrix');
			end
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function set.delta(corr_map,delta)
			
			if (corr_map.delta ~= delta)
				corr_map.delta = delta;
				notify(corr_map,'outdated_corr_matrix');
			end
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function set.filter_type(corr_map,filter_type)
			
			if ~strcmp(corr_map.filter_type,filter_type)
				corr_map.filter_type = filter_type;
				notify(corr_map,'outdated_corr_matrix');
			end
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function H = plot_corr_signal(corr_map,axes_handle,varargin)
			
			% Plot the 1D correlation signal
			H = plot(axes_handle,corr_map.time_vector,corr_map.corr_signal,varargin{:});
			
			% Set axes limits
			% xlim(axes_handle,corr_map.time_limits);
			ylim(axes_handle,[-1 +1]);
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		% 		function H = plot_corr_map(corr_map,axes_handle)
		%
		% 			% The x and y limits, in seconds, of the the correlation map
		% 			x_lim = corr_map.time_limits;
		% 			y_lim = [corr_map.time_offsets(end) corr_map.time_offsets(1)];
		%
		% 			% Raise the colormap resolution (the default value is 64 colors)
		% 			colormap(jet(100));
		%
		% 			% Plot the correlation map
		% 			% Important: For some weird reason, if we specify the parent axes
		% 			% when calling the function imagesc(), the range of the associated
		% 			% colorbar will not be [-1 +1] as specified in the function caxis()
		% 			% below, but will instead be [1 100] (which probably comes from the
		% 			% number of entries in the colormap, as defined above). What we do
		% 			% then is to set the current axes to 'axes_handle' and then call
		% 			% the imagesc() function without specifying the parent axes.
		% 			% H = imagesc(x_lim,y_lim,corr_map.corr_map,'Parent',axes_handle);
		% 			axes(axes_handle); H = imagesc(x_lim,y_lim,corr_map.corr_map);
		% 			set(axes_handle,'YDir','normal');
		%
		% 			% Set the colormap range to [-1 +1]. This means that dark blue
		% 			% and dark red will map to -1 and +1, respectively, and not to
		% 			% the minimum and maximum values in the data
		% 			caxis(axes_handle,[-1 +1]);
		%
		% 		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function [H1,H2] = plot_corr_map(corr_map,axes_handle,colorbar_position)
			
			% Plot the correlation map on the axes with handle 'axes_handle'.
			% If 'colorbar_position' is provided, a colorbar is created at
			% the specified position. The colorbar position is considered
			% to be expressed in the same units as the provided axes.
			% Although Matlab can resize the axes associated with a colorbar
			% in order to accommodate it, this method won't allow that. In
			% other words, the position of the provided axes won't change.
			% Thus, the positions of both the colormap and the colorbar axes
			% should be defined appropriatelly before this method is called.
			
			
			% The x and y limits, in seconds, of the the correlation map
			x_lim = corr_map.time_limits;
			y_lim = [corr_map.time_offsets(end) corr_map.time_offsets(1)];
			
			% Raise the colormap resolution (the default value is 64 colors)
			colormap(jet(100));
			
			% Plot the correlation map
			% Important: For some weird reason, if we specify the parent axes
			% when calling the function imagesc(), the range of the associated
			% colorbar will not be [-1 +1] as specified in the function caxis()
			% below, but will instead be [1 100] (which probably comes from the
			% number of entries in the colormap, as defined above). What we do
			% then is to set the current axes to 'axes_handle' and then call
			% the imagesc() function without specifying the parent axes.
			% H = imagesc(x_lim,y_lim,corr_map.corr_map,'Parent',axes_handle);
			axes(axes_handle); H1 = imagesc(x_lim,y_lim,corr_map.corr_map);
			set(axes_handle,'YDir','normal');
			
			% Set the colormap range to [-1 +1]. This means that dark blue
			% and dark red will map to -1 and +1, respectively, and not to
			% the minimum and maximum values in the data
			caxis(axes_handle,[-1 +1]);
			
			% If a colorbar position has been provided, we go ahead and create
			% the colorbar
			if exist('colorbar_position','var')
				
				% The requirements for the input argument 'colorbar_position'
				requirements(1) = isvector(colorbar_position);
				requirements(2) = isnumeric(colorbar_position);
				requirements(3) = (length(colorbar_position) == 4);
				
				% All requirements must be met
				if ~all(requirements)
					error('The colorbar position must be a four-element numeric vector');
				end
				
				% The units and position of the colormap axes
				axes_units = get(axes_handle,'Units');
				axes_position = get(axes_handle,'Position');
				
				% Create the colorbar at the position specified by the variable
				% 'colorbar_position' and associate it with the colormap axes
				H2 = colorbar('Peer',axes_handle,'Units',axes_units,'Position',colorbar_position);
				
				% Restore the position of the colormap axes. This ensures the
				% colormap axes are at the same position they were before the
				% colorbar was created. This is necessary because the colormap
				% axes may have been resized by the function colorbar() in
				% order to accommodate the colorbar.
				set(axes_handle,'Position',axes_position);
				
			end
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
	end
	
	% --------------------------------------------------------------------------------------------- %
	
	methods (Access = private)
		
		% --------------------------------------------------------------------------------------------- %
		
		function set_signals(corr_map,signal_1,signal_2)
			
			% Both input signals must be 'TimeSignal' objects
			if ~(isa(signal_1,'TimeSignal') & isa(signal_2,'TimeSignal'))
				error('Input signals must be ''TimeSignal'' objects.');
			end
			
			% Input signals must have the same rate
			rate = unique([signal_1.rate signal_2.rate]);
			if isscalar(rate)
				corr_map.rate = rate;
			else
				error('Input signals must have the same rate.');
			end
			
			% Input signals must have the same number of samples
			n_samples = unique([signal_1.n_samples signal_2.n_samples]);
			if isscalar(n_samples)
				corr_map.n_samples = n_samples;
			else
				error('Input signals must have the same number of samples.');
			end
			
			% The input signals
			corr_map.signal_1 = signal_1;
			corr_map.signal_2 = signal_2;
			
			% The names of the input signals
			% corr_map.input_signal_names{1} = signal_1.name;
			% corr_map.input_signal_names{2} = signal_2.name;
			
			% We use signal 1 as the reference signal, and so the sample
			% is taken as the start sample of signal 1
			corr_map.start_sample = signal_1.start_sample;
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function corr_matrix_to_corr_map(corr_map,mode)
			
			% If 'mode' is not provided, it defaults to 'full'
			if (nargin < 2), mode = 'full'; end;
			
			% The valid modes
			valid_modes = {'half' 'full'};
			
			% Check if 'mode' is a valid mode
			if ~ismember(mode,valid_modes)
				error_message = sprintf('Invalid mode: ''%s''. Valid modes are ''half'' and ''full''.',mode);
				error(error_message);
			end
			
			% The lag range (-delta:delta)
			delta = corr_map.delta;
			
			% A matrix containing the diagonals of the correlation matrix
			% (corr_map.corr_matrix) along its columns
			mapa = spdiags(corr_map.corr_matrix,[-delta:delta]);
			
			% Number of rows and columns in the correlation map matrix
			[n_rows n_cols] = size(mapa);
			
			% The middle column
			mid_col = delta+1;
			
			% An auxiliary matrix
			% TODO: FINISH COMMENTING THIS METHOD!!!
			aux = zeros(2*n_rows,n_cols);
			aux(1:2:end) = mapa;
			
			h = [0.5 1 0.5];
			aux = filter(h,1,aux);
			aux(1,:) = [];
			
			mapa = zeros(size(aux));
			
			for m=0:mid_col-1
				mapa(1+m:end-m,mid_col-m) = aux(1:end-2*m,mid_col-m);
				mapa(1+m:end-m,mid_col+m) = aux(1+2*m:end,mid_col+m);
			end
			
			% If 'mode' is 'half', we decimate the correlation map
			% by a factor of 2 along the time dimension
			if isequal(mode,'half'), mapa = mapa(1:2:end,:); end;
			
			% Rotate the correlation map 90 degrees (so that time is along the
			% horizontal direction) and store it in the CorrelationMap object
			corr_map.corr_map = rot90(mapa);
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function update_corr_matrix(corr_map,evento)
			
			% disp('Updating the correlation matrix...');
			
			corr_map.corr_map    = [];
			corr_map.corr_matrix = [];
			corr_map.corr_signal = [];
			
			if (corr_map.delta > 0)
				corr_map.corr_matrix = instantaneous_correlation(corr_map.signal_1.signal,corr_map.signal_2.signal,corr_map.filter_type,corr_map.eta,corr_map.delta);
				corr_map.corr_signal = full(diag(corr_map.corr_matrix));
				% corr_map.corr_map = rot90(spdiags(corr_map.corr_matrix));
				corr_map.corr_matrix_to_corr_map();
			else
				corr_map.corr_signal = instantaneous_correlation(corr_map.signal_1.signal,corr_map.signal_2.signal,corr_map.filter_type,corr_map.eta,0);
			end
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
	end
	
	% --------------------------------------------------------------------------------------------- %
	
end
