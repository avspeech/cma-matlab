function [figure_handle,axes_handles] = create_panel_layout(panel_widths,panel_heights,panel_spacing,margins,figure_handle)

% [FIGURE_HANDLE,AXES_HANDLES] = CREATE_PANEL_LAYOUT(PANEL_WIDTHS,PANEL_HEIGHTS,PANEL_SPACING,MARGINS,FIGURE_HANDLE)
% creates a figure with multiple panels (axes) distributed in a grid according to the provided layout.
% PANEL_WIDTHS and PANEL_HEIGHTS are vectors specifying the widths and heights of the panels, respectively.
% The number of elements in PANEL_WIDTHS and PANEL_HEIGHTS determine the number columns and rows in the panel
% grid, respectively. Thus, PANEL_WIDTHS(I) determines the width of the panels in the I-th column of the grid,
% whereas PANEL_HEIGHTS(J) determines the height of the panels in the J-th row of the grid. PANEL_SPACING is
% a two element vector of the form [horz_spacing vert_spacing], where 'horz_spacing' and 'vert_spacing' are
% the horizontal and vertical spacing between the panels, respectively. MARGINS is a four element vector of
% the form [left_margin right_margin top_margin bottom_margin] with the size of the figure margins. Lastly,
% FIGURE_HANDLE is a positive integer which will be used as the figure handle. If a figure with the provided
% handle already exists, an error occurs. All provided sizes must be in centimeters.
% 
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2009-2014 Adriano Vilela Barbosa


% All input variables in centimeters!!!
% panel_spacing = [horz_spacing vert_spacing]
% margins       = [left right top bottom]

% The number of elements in 'panel_widths' determines the number of panel columns
% The number of elements in 'panel_heights' determines the number of panel rows

% Non provided input arguments default to an empty matrix to indicate
% that a new figure must be created
if (nargin < 5), figure_handle = []; end;
if (nargin < 4), margins       = []; end;
if (nargin < 3), panel_spacing = []; end;

% The number of rows and columns in the "panel matrix"
n_rows = length(panel_heights);
n_cols = length(panel_widths);

% If the margins are not provided, they default to 2 cm
if isempty(margins), margins = [2 2 2 2]; end;

% If the panel spacing is not provided, it defaults to 1 cm
if isempty(panel_spacing), panel_spacing = [1 1]; end;

% The figure size, computed from the panel sizes, inter-panel spacings, and margins.
% The figure size is a two-element vector of the form [figure_width figure_height]
figure_size(1) = sum(margins(1:2)) + sum(panel_widths) + (n_cols-1)*panel_spacing(1);
figure_size(2) = sum(margins(3:4)) + sum(panel_heights) + (n_rows-1)*panel_spacing(2);

% Create figure
figure_handle = create_figure(figure_size,'Centimeters',figure_handle);

% Reorganize the panel heights so that they are listed from bottom to top
panel_heights = flipud(panel_heights(:));

% Panels are number from left to right, bottom to top
for m=1:n_rows

	% The vertical position of the lower edge of the current panel
	vert_position = margins(4) + (m-1)*panel_spacing(2) + sum(panel_heights(1:m-1));

	for n=1:n_cols

		% The horizontal position of the left edge of the current panel
		horz_position = margins(1) + (n-1)*panel_spacing(1) + sum(panel_widths(1:n-1));

		% The position of lower left edge of the current panel
		panel_corner = [horz_position vert_position];

		% The position of the current panel (in Matlab notation)
		panel_position = [panel_corner panel_widths(n) panel_heights(m)];

		% Create the current panel at the specified position. It's very
		% important that we set the axes units before setting its size,
		% size the latter depends on the former
		axes_handles(m,n) = axes('Units','Centimeters','Position',panel_position);

	end

end

% Reorganize the axes handles so that the first handle corresponds
% to the top panel and the last handle to the bottom panel
% handles.axes = flipud(handles.axes);
axes_handles = flipud(axes_handles);

% Set both the figure's and the axes' units to 'normalized'
set(figure_handle,'Units','Normalized');
set(axes_handles,'Units','Normalized');

% ------------------------------------------------------------------------------------------------ %
