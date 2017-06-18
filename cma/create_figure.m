function figure_handle = create_figure(figure_size,figure_units,figure_handle)

% HANDLE = CREATE_FIGURE(SIZE,UNITS,HANDLE) creates a figure with handle
% HANDLE, of size SIZE, measured in units UNITS. SIZE = [width height]
% is a two element vector containing the figure width and height, in
% this order. UNITS can be any valid Matlab unit (see the Matlab help
% on figure properties for details). HANDLE is a positive integer which
% will be used as the figure handle. If a figure with the provided handle
% already exists, an error occurs. If SIZE is not provided, the default
% size [16 12] is used. If UNITS is not provided, centimeters are used.
% If HANDLE is not provided, the figure handle is automatically chosen
% by Matlab. The returned argument HANDLE contains the handle to the
% created figure.
%
% This function automatically centers the created figure on the screen,
% making it unnecessary to specify the lower left position of the figure,
% which is the case when calling Matlab's function figure() directly.
% This is specially useful when dealing with un-normalized units.
%
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2008-2014 Adriano Vilela Barbosa


% Default values
default_units	= 'centimeters';
default_size	= [16 12];
default_handle = [];

% Non provided properties are assigned default values
if (nargin < 3), figure_handle = default_handle; end;
if (nargin < 2), figure_units = default_units; end;
if (nargin < 1), figure_size = default_size; end;

% All figure handles
handle_list = findobj(0,'Type','figure','-depth',1);

% Check if a figure with the provided handle already exists
if ismember(figure_handle,handle_list)
	error(sprintf('A figure with the provided handle (%g) already exists.',figure_handle));
end

% Argument FIGURE_SIZE must be a two element vector
if ~all([isvector(figure_size) (length(figure_size) == 2)])
	error('Argument FIGURE_SIZE must be a two element vector.');
end

% Ensures 'figure_size' is a row vector
figure_size = figure_size(:)';

% The screen size, expressed in the same units as the figure
aux = get(0,'Units');
set(0,'Units',figure_units);
screen_size = get(0,'ScreenSize');
screen_size = screen_size(3:4);
set(0,'Units',aux);

% The position where the new figure will be placed
margins = (screen_size-figure_size)/2;
position = [margins figure_size];

% Create new figure
if isempty(figure_handle)
	figure_handle = figure();
else
	figure_handle = figure(figure_handle);
end

% Set the position of the new figure
set(figure_handle,'Units',figure_units,'Position',position);

%-------------------------------------------------------------------------------------------------%
