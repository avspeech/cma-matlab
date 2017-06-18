function align_labels(figure_handle)

% ALIGN_LABELS(FIGURE_HANDLE) aligns the ylabels of the various axes
% of the figure FIGURE_HANDLE. If the figure handle is not provided,
% the current figure is used.
% 
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2003-2014 Adriano Vilela Barbosa


% If the figure handle is not provided, we use gcf
if (nargin<1), figure_handle = gcf; end;

% The axes handles
% axes_handles = sort(findobj(figure_handle,'Type','axes'));
axes_handles = sort(findobj(figure_handle,'Type','axes','Tag',''));

% The ylabel handles
for m=1:length(axes_handles)
	y_labels(m) = get(axes_handles(m),'YLabel');
end

% The ylabel positions
for m=1:length(y_labels)

	% The default units for the label positions is 'Data'. We cannot
	% ensure that these units are the same across the different axes
	% of the figure. Thus, before proceeding, we express all label
	% positions in a common units system
	% set(y_labels(m),'Units','Normalized');
	set(y_labels(m),'Units','Centimeters');

	% The ylabel positions, in normalized units
	ylabel_positions(m,:) = get(y_labels(m),'Position');

end

% The horizontal position of the leftmost ylabel
ylabel_positions(:,1) = min(ylabel_positions(:,1));

% Moves all ylabels to the same horizontal position
for m=1:length(y_labels)
	set(y_labels(m),'Position',ylabel_positions(m,:));
end

%-------------------------------------------------------------------------------------------------%
