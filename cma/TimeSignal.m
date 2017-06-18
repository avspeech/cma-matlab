classdef TimeSignal < handle

	% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
	% Copyright 2009-2014 Adriano Vilela Barbosa
	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
	% Public properties
	properties
		signal       = [];
		name         = '';
		start_sample = 1;
	end
	
	% Private, non-dependent properties
	properties(SetAccess = private)
		rate = 1;
	end
	
	% Private, dependent properties
	% Update on 17-feb-2012: we have renamed the property 'time_limits' to
	% 'time_range'. However, we have kept the old property for the sake of
	% backward compatibility
	properties(Dependent = true, SetAccess = private)
		n_samples   = [];
		time_vector = [];
		time_limits = [];
		time_range  = [];
	end
	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
	methods
		
		% TODO: WE MAY WANT TO IMPLEMENT A METHOD CALLED 'EXTRACT SEGMENT'
		% WHICH RETURNS A NEW TIMESIGNAL OBJECT CONTAINING ONLY THE SEGMENT
		% OF INTEREST. THIS METHOD COULD TAKE AN OPTIONAL ARGUMENT TELLING
		% US IF THE TIME BOUNDARIES OF THE NEW OBJECT SHOULD BE RELATIVE
		% TO ITS PARENT OBJECT OR NOT.
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_signal = TimeSignal(signal,rate,start_sample)
			
			% This is the constructor function
			
			% Handle the zero argument case that can happen, for example,
			% when a subclass of the current class is instantiated
			if (nargin == 0), return; end;
			
			% Initialize the signal
			time_signal.signal = signal;
			
			% Fill in the remaining properties, if they have been provided
			
			
			% Important: the property 'rate' has a setter that updates the property
			% 'start_sample' everytime the signal rate changes. So, whenever assigning
			% values to both 'rate' and 'start_sample', they have to be updated in the
			% correct order ('rate' and then 'start_sample'). If we set 'start_sample'
			% first and then 'rate', the setter for the latter will overwrite the value
			% of the former. This problem happens when loading TimeSignal objects from
			% Matlab mat files. I think the easiest way of fixing is just to get rid
			% of the 'start_sample' property and use a 'start_time' property instead.
			% Update on 17-feb-2012: We got rid of the setter for the property 'rate'.
			% The idea behind using this setter in the first place was to be able to
			% resample signals by simply setting their 'rate' property to a different
			% value. However, this has turned out to be quite problematic. If we need
			% it, we can always add a resampling method later which will update both
			% the 'rate' and 'start_sample' properties after resampling the signal.
			
			% if (nargin > 1), time_signal.rate = rate; end;
			if (nargin > 2), time_signal.start_sample = start_sample; end;
			if (nargin > 1), time_signal.rate = rate; end;
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function n_samples = get.n_samples(time_signal)
			n_samples = length(time_signal.signal);
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_vector = get.time_vector(time_signal)
			
			% The time vector
			time_vector = (time_signal.start_sample:time_signal.start_sample+time_signal.n_samples-1)/time_signal.rate;
			time_vector = time_vector(:);
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_limits = get.time_limits(time_signal)
			time_vector = time_signal.time_vector;
			time_limits = [time_vector(1) time_vector(end)];
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function time_range = get.time_range(time_signal)
			time_vector = time_signal.time_vector;
			time_range = [time_vector(1) time_vector(end)];
		end
		
		% --------------------------------------------------------------------------------------------- %
		
		% We have removed the setter for the property 'rate'. See comments on
		% the constructor method for the reasons we decided to do that.
		
		% 		function set.rate(time_signal,new_rate)
		%
		% 			% TODO: WE HAVE CHANGED THE PROPERTY 'RATE' FROM PUBLIC TO PRIVATE.
		% 			% NOW WE HAVE TO CONVERT THIS FUNCTION INTO A RESAMPLING METHOD.
		%
		% 			% This function changes the property 'rate' and then
		% 			% updates the property 'start_sample' accordingly
		%
		% 			% The current signal rate
		% 			old_rate = time_signal.rate;
		%
		% 			% The ratio between the new and the old rates
		% 			q = new_rate/old_rate;
		%
		% 			% Change the signal rate to its new value
		% 			time_signal.rate = new_rate;
		%
		% 			% Express the start sample at the new rate. If the resampling
		% 			% operation results in a non-integer start sample, we round
		% 			% it towards the nearest integer
		% 			time_signal.start_sample = round((time_signal.start_sample-1)*q+1);
		%
		% 		end
		
		% --------------------------------------------------------------------------------------------- %
		
		function H = plot(time_signal,axes_handle,varargin)
			
			% Non-provided input arguments are assigned an empty matrix
			% to indicate that they should get their default values
			if (nargin < 2), axes_handle = []; end;
			
			% If no axes where to plot the signal has been provided,
			% we plot it to the current axes
			if isempty(axes_handle), axes_handle = gca; end;
			
			% Plot the time signal
			% H = plot(axes_handle,time_signal.time_vector,time_signal.signal);
			H = plot(axes_handle,time_signal.time_vector,time_signal.signal,varargin{:});
			
			% Labeling
			% xlabel('Time (s)');
			
		end
		
		% --------------------------------------------------------------------------------------------- %
		
	end
	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
% 	methods (Access = protected)
% 		
% 		% --------------------------------------------------------------------------------------------- %
% 		
% 		function signal_mat = enframe_signal(time_signal,varargin)
% 			
% 			% signal_mat = enframe_speech(signal_vec,frame_length,frame_shift,padding);
% 			signal_mat = enframe_speech(time_signal.signal,varargin{:});
% 			
% 		end
% 		
% 		% --------------------------------------------------------------------------------------------- %
% 		
% 	end
	
	% --------------------------------------------------------------------------------------------- %
	% --------------------------------------------------------------------------------------------- %
	
end


% --------------------------------------------------------------------------------------------- %
