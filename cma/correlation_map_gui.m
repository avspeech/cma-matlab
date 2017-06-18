function [handles,corr_maps] = correlation_map_gui(signal_1,signal_2,audio_signals,filter_type,eta,delta,plot_1d,interactive,plot_inputs_together)

% [HANDLES,CORR_MAPS] = CORRELATION_MAP_GUI(SIGNAL_1,SIGNAL_2,AUDIO_SIGNALS,FILTER_TYPE,ETA,DELTA,PLOT_1D,INTERACTIVE,PLOT_INPUTS_TOGETHER)
% computes and plots the correlation map between signals SIGNAL_1 and SIGNAL_2.
% 
% 
% INPUTS:
% 
% signal_1, signal_2:	the two signals to be correlated. These are the only mandatory input arguments.
% audio_signals:			signals provided here will be plotted along with the two input signals, in separate panels. This is useful when
%								correlating motion signals and we want to visualize the corresponding audio signals. This is primarily meant for
%								audio signals, but any signals can be used.
% filter_type:				either 'ss' (single-sided, uni-directional) or 'ds' (double-sided, bi-directional)
% eta:						the exponential decay, which controls the filter sensitivity. It must be a small positive number (between zero
%								and one). Important: 'eta' can be a vector. In that case, one correlation panel is plotted for each value of eta.
% delta:						the lag range of the correlation map, in samples (range = [-delta +delta])
% plot_1d:					if true, the 1D correlation signal (which is simply the 2D correlatin map at lag zero) will be plotted.
% interactive:				if true, widgets for controlling the filter type and the values of eta and delta will be available
% plot_inputs_together:	if true, the two input signals will be plotted in the same panel
% 
% OUTPUTS:
% 
% handles:		a structure with handles to the various panels. This is useful for tweaking the plots later, if necessary.
% corr_maps:	a structure with the correlation maps
% 
% ---
% 
% Note: only the first two input arguments are mandatory. All other arguments, if not provided, or provided as empty matrices,
% will be assigned their default values
% 
% 
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2006-2014 Adriano Vilela Barbosa


% Non provided input arguments are assigned an empty matrix, which
% means they will be assigned their default values somewhere else
if (nargin<9), plot_inputs_together = []; end;
if (nargin<8), interactive          = []; end;
if (nargin<7), plot_1d              = []; end;
if (nargin<6), delta                = []; end;
if (nargin<5), eta                  = []; end;
if (nargin<4), filter_type          = []; end;
if (nargin<3), audio_signals        = []; end;

% If 'plot_inputs_together' is not provided, the two input signals
% are plotted in separate panels.
if isempty(interactive), plot_inputs_together = false; end;

% If 'interactive' is not provided, it defaults to 'false'. Important:
% do not change this default value, at the risk of breaking old code.
% At some point, this function was non-interactive only, and the input
% argument 'interactive' didn't exist. In order to keep backward com-
% patibility with code using that old interface, the flag 'interactive'
% must default to false when not provided.
if isempty(interactive), interactive = false; end;

% If 'plot_1d' is not provided, it defaults to 'false'
if isempty(plot_1d), plot_1d = false; end;

% If 'delta' is not provided, it defaults to zero
if isempty(delta), delta = 0; end;

% If 'eta' is not provided, it defaults to 0.1
if isempty(eta), eta = 0.1; end;

% If 'filter_type' is not provided, it defaults to 'ds'
if isempty(filter_type), filter_type = 'ds'; end;

% Interactive plots require a scalar 'eta'
if (interactive & ~isscalar(eta))
	error('Interactive plots require a scalar ''eta''.');
end

% Compute the instantaneous correlation between the two signals
for m=1:length(eta)
	corr_maps(m) = CorrelationMap(signal_1,signal_2,filter_type,eta(m),delta);
end

% Create the GUI and plot the signals
% handles = create_gui(audio_signal,corr_maps,plot_1d);
% handles = create_gui(audio_signals,corr_maps,plot_1d,interactive);
handles = create_gui(audio_signals,corr_maps,plot_inputs_together,plot_1d,interactive);

% Organizes the data into a structure that will be shared among functions
dados.audio_signals = audio_signals;
dados.corr_maps     = corr_maps;
dados.handles       = handles;

% Stores the shared structure 'dados'
setappdata(handles.figure,'dados',dados);

% --------------------------------------------------------------------------------------------- %

% function handles = create_gui(audio_signal,corr_maps,plot_1d)
% function handles = create_gui(audio_signals,corr_maps,plot_1d,interactive)
function handles = create_gui(audio_signals,corr_maps,plot_inputs_together,plot_1d,interactive)


% Create the graphical user interface and plot the signals

% Possibilites:
%
% delta = 0			plot_1d = false			=> plot_1d = true; plot_2d = false
% delta = 0			plot_1d = true			=> plot_1d = true; plot_2d = false
% delta ~= 0		plot_1d = false			=> plot_1d = false; plot_2d = true
% delta ~= 0		plot_1d = true			=> plot_1d = true;  plot_2d = true



% Is there an audio signal?
audio_present = ~isempty(audio_signals);

% The number of audio signals
n_audio_signals = length(audio_signals);

% The number of eta values
n_etas   = length(corr_maps);

% The offset range 'delta'. Please note that all CorrelationMap objects
% have the same value of 'delta'
delta = corr_maps(1).delta;

% The two input signals
signal_1 = corr_maps(1).signal_1;
signal_2 = corr_maps(1).signal_2;

% The boolean 'plot_2d' is initialized to true
plot_2d = true;

% If 'delta' is zero, 'plot_2d' is set to false (as there is
% no 2D correlation map to plot) and 'plot_1d' is set to true
% (as this is the only correlation signal available for plotting)
if (delta == 0), plot_1d = true; plot_2d = false; end;

% The font size to be used for labeling axes
label_font_size = 1.1*get(0,'DefaultAxesFontSize');

% Panel widths, in centimeters
panel_widths.main_panels = 16;
panel_widths.colorbars   = 0.6;

% Panel ratio (the ratio between the heights of the signal panels
% and the correlation map panels)
panel_ratio = 2.2;

% Panel heights, in centimeters
panel_heights.signals   = 2.4;
panel_heights.corr_maps = panel_ratio*panel_heights.signals;

% If no correlation maps are being plotted, we set the dimensions
% of the axes associated with them to empty matrices.
if ~plot_2d
	panel_widths.colorbars  = [];
	panel_heights.corr_maps = [];
end

% Horizontal and vertical spacing between panels, in centimeters
panel_spacing = [0.5 0.4];

% The figure margins, in centimeters (left, right, top, bottom)
% margins = [2 1 1 2];
margins = [2.0 1.5 1.0 2.2];


% TODO: IT WOULD BE VERY USEFUL TO BE ABLE TO PASS THE PANEL SIZES
% TO THE correlation_map_gui() FUNCTION. I FOUND MYSELF EDITING THIS
% FILE ALL THE TIME IN ORDER TO GET THE PANELS THE EXACT SIZE I WANT.
% THE IMPORTANT VARIABLE HERE IS 'panel_heights.signals'. THE HEIGHTS
% OF THE CORRELATION MAP PANELS ARE COMPUTED FROM THE HEIGHTS OF THE
% SIGNAL PANELS USING THE VARIABLE 'panel_ratio'.


% The position, in centimeters, of the widgets panel (the panel used
% as a container for the GUI widgets)
widgets_panel_position = [margins(1) 1 panel_widths.main_panels 3];

% Extra margins to be added to the figure margins. These extra margins
% are initialized to zero but are increased in case of interactive plots,
% due to the extra space required for the widgets panel
extra_margins = [0 0 0 0];

% If the plot is interactive, we add extra space to the bottom
% margin of the figure in order to place the GUI widgets there
if interactive
	extra_margins(end) = sum(widgets_panel_position([2 4]));
end

% The number of panels for the correlation maps
n_panels.corr_maps = n_etas*length(panel_heights.corr_maps);

% The number of panels for the various 1D signals (audio, input signals,
% 1D correlation signals)
% n_panels.signals = 2;
if plot_inputs_together, n_panels.signals = 1; else, n_panels.signals = 2; end;
% if audio_present, n_panels.signals = n_panels.signals+1; end;
if audio_present, n_panels.signals = n_panels.signals+n_audio_signals; end;
if plot_1d, n_panels.signals = n_panels.signals+n_etas; end;

% Panel widths and heights, as expected by the function create_panel_layout()
panel_widths = [panel_widths.main_panels panel_widths.colorbars];
panel_heights = ...
	[repmat(panel_heights.signals,1,n_panels.signals) ...
	repmat(panel_heights.corr_maps,1,n_panels.corr_maps)];

% Create a figure with all the necessary axes
% [figure_handle,axes_handles] = create_panel_layout(panel_widths,panel_heights,panel_spacing,margins);
[figure_handle,axes_handles] = create_panel_layout(panel_widths,panel_heights,panel_spacing,margins+extra_margins);
set(axes_handles,'Box','on','XTick',[],'YTick',[]);

% Initialize the structure 'handles'
handles.figure        = figure_handle;
handles.audio_signals = [];
handles.signal_1      = [];
handles.signal_2      = [];
handles.corr_signals  = [];			% 1D correlations
handles.corr_maps     = [];			% 2D correlations
handles.colorbars     = [];

% The handle for the correlation map image. It will be
% initialized by the function plot_correlation_map()
handles.corr_map_images = [];

% The handles for the marker on the correlation map and for the vertical
% lines in the input signals. These handles are updated by the function
% update_cursor_position()
handles.marker = [];
handles.line_1 = [];
handles.line_2 = [];

% Initialize the panel counter. Panels are numbered from top to bottom.
counter = 1;

% If an audio signal is present, it is associated with the first
% panel and plot the signal
if audio_present
	
	% Loop over audio signals
	for m=1:n_audio_signals
		
		% Assign panel number 'counter' to the audio signal
		handles.audio_signals(m) = axes_handles(counter,1);
		counter = counter+1;
		
		% Plot the audio signal and do some labeling
		plot(handles.audio_signals(m),audio_signals(m).time_vector,audio_signals(m).signal,'k');
		set(handles.audio_signals(m),'XTickLabel',[]);
		ylabel(handles.audio_signals(m),audio_signals(m).name,'FontSize',label_font_size,'Interpreter','None');
		
	end
	
end

% The line colors and widths, as well as the signal names, depend on
% whether the two signals are being plotted together or separately
if plot_inputs_together
	axes_indices = [counter counter];
	line_colors  = {'r' 'k'};
	line_widths  = [1.5 0.5];
	signal_names = sprintf('%s, %s',signal_1.name,signal_2.name);
	signal_names = {signal_names signal_names};
else
	axes_indices = [counter counter+1];
	line_colors  = {'k' 'k'};
	line_widths  = [0.5 0.5];
	signal_names = {signal_1.name signal_2.name};
end

% Associate the next panels with the two input signals
handles.signal_1 = axes_handles(axes_indices(1),1);
handles.signal_2 = axes_handles(axes_indices(2),1);
counter = axes_indices(end)+1;

% Plot input signal 1
H = signal_1.plot(handles.signal_1);
hold(handles.signal_1,'on');
set(H,'Color',line_colors{1},'LineWidth',line_widths(1));
set(handles.signal_1,'XTickLabel',[]);
ylabel(handles.signal_1,signal_names{1},'FontSize',label_font_size,'Interpreter','None');

% Plot input signal 2
H = signal_2.plot(handles.signal_2);
hold(handles.signal_2,'on');
set(H,'Color',line_colors{2},'LineWidth',line_widths(2));
set(handles.signal_2,'XTickLabel',[]);
ylabel(handles.signal_2,signal_names{2},'FontSize',label_font_size,'Interpreter','None');

% Initialize the vertical lines
handles.line_1 = line([0 0],[0 0],'Parent',handles.signal_1,'Color',[0 0 0]);
handles.line_2 = line([0 0],[0 0],'Parent',handles.signal_2,'Color',[0 0 0]);

% If there are 1D correlation signals, associate them with the next
% 'n_etas' panels and plot the signals
if plot_1d
	
	handles.corr_signals = axes_handles(counter:counter+n_etas-1,1);
	counter = counter+n_etas;
	
	for m=1:n_etas
		H = corr_maps(m).plot_corr_signal(handles.corr_signals(m));
		set(H,'Color',[0 0 0]);
		% xlim(handles.corr_signals(m),corr_maps(m).time_limits);
		set(handles.corr_signals(m),'XTickLabel',[]);
		x_lim = xlim(handles.corr_signals(m));
		% axes(handles.corr_signals(m)); line(x_lim,[0 0],'Color',[0 0 0],'LineStyle',':');
		axes(handles.corr_signals(m)); line(x_lim,[0 0],'Color',[0 0 0],'LineStyle',':','HandleVisibility','off');
		ylabel(handles.corr_signals(m),'1D Corr.','FontSize',label_font_size,'Interpreter','None');
	end
	
end

% If there are 2D correlation signals (corr maps), associate them with
% the next 'n_etas' panels and plot the signals
if plot_2d
	
	handles.corr_maps = axes_handles(counter:counter+n_etas-1,1);
	handles.colorbars = axes_handles(counter:counter+n_etas-1,2);
	
	% Delete the unnecessary colorbar axes (these axes were created
	% with the sole purpose of aligning the remaining axes)
	aux_handles = setdiff(axes_handles(:,2),handles.colorbars);
	delete(aux_handles);
	
	for m=1:n_etas
		
		% Plot the correlation map
		colorbar_position = get(handles.colorbars(m),'Position');
		delete(handles.colorbars(m));
		[handles.corr_map_images(m),handles.colorbars(m)] = corr_maps(m).plot_corr_map(handles.corr_maps(m),colorbar_position);
		
		% Do some labeling
		set(handles.corr_maps(m),'XTickLabel',[]);
		% aux = sprintf('Offset (s), eta = %1.2f',corr_maps(m).eta);
		aux = 'Offset (s)';
		if (n_etas>1), aux = sprintf('%s, eta = %1.2f',aux,corr_maps(m).eta); end;
		ylabel(handles.corr_maps(m),aux,'FontSize',label_font_size);
		
		% Initialize the correlation map marker
		handles.marker(m) = text(0,0,'+','Parent',handles.corr_maps(m),'Visible','Off','HitTest','Off','HorizontalAlignment','Center','FontSize',18);
		
		% Assign a callback function to the correlation map image
		set(handles.corr_map_images(m),'ButtonDownFcn',{@update_cursor_position,handles.figure});
		
	end
	
end

% The handle for the bottom-most panel
if plot_2d
	last_axes_handle = handles.corr_maps(end);
else
	last_axes_handle = handles.corr_signals(end);
end

% Some x-labeling for the bottom-most panel
set(last_axes_handle,'XTickLabelMode','Auto');
xlabel(last_axes_handle,'Time (s)','FontSize',label_font_size,'Interpreter','None');

% Align ylabels across all axes objects
align_labels(handles.figure);

% Link all axes in the x direction
% linkaxes([handles.audio_signal; handles.signal_1; handles.signal_2; handles.corr_signals(:); handles.corr_maps(:)],'x');
% linkaxes([handles.corr_maps(:); handles.corr_signals(:); handles.audio_signals(:); handles.signal_1; handles.signal_2],'x');
linkaxes(unique([handles.corr_maps(:); handles.corr_signals(:); handles.audio_signals(:); handles.signal_1; handles.signal_2]),'x');

% Ensure the figure toolbar is shown. The toolbar is very useful
% for zooming in/out the plots
set(handles.figure,'Toolbar','figure');

% If the plot isn't interactive, there's nothing left to do and we then
% return to the calling function. Otherwise, we go ahead and create the
% GUI widgets
if ~interactive, return; end;


% ----- Create widgets for interaction with the plots ----- %


% We define widgets for controlling the following properties:
% eta         -> slider and edit text box
% delta       -> edit text box
% filter_type -> popup menu


handles.gui.widgets_panel = uipanel(handles.figure,...
	'Units','Centimeters','Position',widgets_panel_position) ;

set(handles.gui.widgets_panel,'Units','Normalized');


% Positions for the static text widgets
positions = [0.05 0.70 0.20 0.20; 0.50 0.70 0.20 0.20; 0.75 0.70 0.20 0.20];
positions = align(positions,'Distribute','HorizontalAlignment');


handles.gui.static_text(1) = uicontrol(handles.gui.widgets_panel,'Style','Text',...
	'Units','Normalized',...
	'Position',positions(1,:),...
	'BackgroundColor',[0.8 0.8 0.8],...
	'String','Value of eta',...
	'HorizontalAlignment', 'Center');

handles.gui.static_text(2) = uicontrol(handles.gui.widgets_panel,'Style','Text',...
	'Units','Normalized',...
	'Position',positions(2,:),...
	'BackgroundColor',[0.8 0.8 0.8],...
	'String','Offset range (sec.)',...
	'HorizontalAlignment', 'Center');

handles.gui.static_text(3) = uicontrol(handles.gui.widgets_panel,'Style','Text',...
	'Units','Normalized',...
	'Position',positions(3,:),...
	'BackgroundColor',[0.8 0.8 0.8],...
	'String','Filter type',...
	'HorizontalAlignment', 'Center');






% Get the value of 'eta' from the CorrelationMap object and format it to
% a string so that we can use it to initialize the corresponding widget
[eta,eta_string] = format_string(corr_maps(1).eta);


indice = 1;
aux = positions(indice,:);
aux(2) = 0.4; aux(3) = 0.1;
aux = [positions(indice,:); aux];
aux = align(aux,'Center','HorizontalAlignment');
aux = aux(2,:);

% Create an editable text box for entereing the value of eta
handles.gui.edit_text(1) = uicontrol(handles.gui.widgets_panel,'Style','Edit',...
	'Units','Normalized',...
	'Position',aux,...
	'BackgroundColor',[0.8 0.8 0.8],...
	'Value',eta,...
	'String',eta_string,...
	'TooltipString','Set the eta value',...
	'Callback',{@on_changed_eta,handles.figure});




% Get the value of 'delta' (in seconds) from the CorrelationMap object
% and format it to a string so that we can use it to initialize the
% corresponding widget
[delta_sec,delta_sec_string] = format_string(corr_maps(1).delta/corr_maps(1).rate);



indice = 2;
aux = positions(indice,:);
aux(2) = 0.4; aux(3) = 0.1;
aux = [positions(indice,:); aux];
aux = align(aux,'Center','HorizontalAlignment');
aux = aux(2,:);

% Create an editable text box for entereing the value of eta
handles.gui.edit_text(2) = uicontrol(handles.gui.widgets_panel,'Style','Edit',...
	'Units','Normalized',...
	'Position',aux,...
	'BackgroundColor',[0.8 0.8 0.8],...
	'Value',delta_sec,...
	'String',delta_sec_string,...
	'TooltipString','Set the delta value',...
	'Callback',{@on_changed_delta,handles.figure});





% The filter type strings. The first cell array contains the strings used
% in the instantaneous correlation functions ('ss' and 'ds'), whereas the
% second cell array contains the strings used in the popup menu widget.
% It is crucial that, for any given element 'n', filter_type_strings_1{n}
% and filter_type_strings_2{n} correspond to the same filter type.
filter_type_strings_1 = {'ss' 'ds'};
filter_type_strings_2 = {'Uni-directional' 'Bi-directional'};

% The value of the current filter type. This will be used to set the
% property 'Value' in the filter type popup menu
[aux,filter_type_value] = ismember(corr_maps.filter_type,filter_type_strings_1);

indice = 3;
aux = positions(indice,:);
aux(2) = 0.4;
% aux(3) = 0.1;
aux = [positions(indice,:); aux];
aux = align(aux,'Center','HorizontalAlignment');
aux = aux(2,:);

% Create a popup menu for choosing the filter type
handles.gui.popup_menu = uicontrol(handles.gui.widgets_panel,'Style','popupmenu',...
	'Units','Normalized',...
	'Position',aux,...
	'BackgroundColor',[0.8 0.8 0.8],...
	'String',filter_type_strings_2,...
	'Value',filter_type_value,...
	'HorizontalAlignment', 'Left',...
	'Callback',{@on_changed_filter_type,handles.figure});




indice = 1;
aux = positions(indice,:);
aux(2) = 0.1;
% aux(3) = 0.1;
aux = [positions(indice,:); aux];
aux = align(aux,'Center','HorizontalAlignment');
aux = aux(2,:);

% Create a slider for controlling the value of eta
handles.gui.slider = uicontrol(handles.gui.widgets_panel,'Style','Slider',...
	'Units','Normalized',...
	'Value',eta,...
	'Position',aux,...
	'TooltipString','Set the eta value',...
	'Callback',{@on_changed_eta,handles.figure});


% --------------------------------------------------------------------------------------------- %

function on_changed_eta(h_object,event_data,figure_handle)

% Load the shared data structure
dados = getappdata(figure_handle,'dados');
handles = dados.handles;
corr_maps = dados.corr_maps;

% Get the new value of eta
switch h_object
	case handles.gui.slider
		% Get the new value of eta from the Slider object
		eta = get(handles.gui.slider,'Value');
	case handles.gui.edit_text(1)
		% Get the new value of eta from the EditText object
		eta = str2num(get(handles.gui.edit_text(1),'String'));
end

% Format string for 'eta'
[eta,eta_string] = format_string(eta);

% Update the CorrelationMap object with the new value of 'eta'
corr_maps.eta = eta;

% Update the Slider object with the new value of 'eta'
set(handles.gui.slider,'Value',eta);
set(handles.gui.slider,'String',eta_string);

% Update the EditText object with the new value of 'eta'
set(handles.gui.edit_text(1),'Value',eta);
set(handles.gui.edit_text(1),'String',eta_string);

% Update the panels that depend on 'eta'
handles = update_corr_plots(handles,corr_maps);

% Stores the shared structure 'dados'
dados.handles = handles;
dados.corr_maps = corr_maps;
setappdata(handles.figure,'dados',dados);

% --------------------------------------------------------------------------------------------- %

function on_changed_delta(h_object,event_data,figure_handle)

% Load the shared data structure
dados = getappdata(figure_handle,'dados');
handles = dados.handles;
corr_maps = dados.corr_maps;

% Get the new value of 'delta' from the EditText object
delta = str2num(get(handles.gui.edit_text(2),'String'));

% Format string for 'delta'
[aux,delta_string] = format_string(delta);

% Update the CorrelationMap object with the new value of 'delta'
corr_maps.delta = ceil(delta*corr_maps.rate);

% Update the EditText object with the new value of 'delta'
set(handles.gui.edit_text(2),'Value',delta);
set(handles.gui.edit_text(2),'String',delta_string);

% Update the panels that depend on 'delta'
handles = update_corr_plots(handles,corr_maps);

% Stores the shared structure 'dados'
dados.handles = handles;
dados.corr_maps = corr_maps;
setappdata(handles.figure,'dados',dados);

% --------------------------------------------------------------------------------------------- %

function on_changed_filter_type(h_object,event_data,figure_handle)

% Load the shared data structure
dados = getappdata(figure_handle,'dados');
handles = dados.handles;
corr_maps = dados.corr_maps;

% The filter type strings. The first cell array contains the strings used
% in the instantaneous correlation functions ('ss' and 'ds'), whereas the
% second cell array contains the strings used in the popup menu widget.
% It is crucial that, for any given element 'n', filter_type_strings_1{n}
% and filter_type_strings_2{n} correspond to the same filter type.
filter_type_strings_1 = {'ss' 'ds'};
filter_type_strings_2 = {'Uni-directional' 'Bi-directional'};

% Get the new value of 'delta' from the EditText object
% delta = str2num(get(handles.gui.edit_text(2),'String'));
filter_string = get(handles.gui.popup_menu(1),'String');
filter_value  = get(handles.gui.popup_menu(1),'Value');

% Update the CorrelationMap object with the new value of 'filter_type'
switch filter_string{filter_value}
	case 'Uni-directional'
		corr_maps.filter_type = 'ss';
	case 'Bi-directional'
		corr_maps.filter_type = 'ds';
end

% % Format string for 'delta'
% [aux,delta_string] = format_string(delta);
%
% % Update the structure 'dados' with the new value of 'delta'
% corr_maps.delta = ceil(delta*corr_maps.rate);
%
% % Update the EditText object with the new value of 'delta'
% set(handles.gui.edit_text(2),'Value',delta);
% set(handles.gui.edit_text(2),'String',delta_string);

% Update the panels that depend on 'delta'
handles = update_corr_plots(handles,corr_maps);

% Stores the shared structure 'dados'
dados.handles = handles;
dados.corr_maps = corr_maps;
setappdata(handles.figure,'dados',dados);

% --------------------------------------------------------------------------------------------- %
%
% function handles = update_corr_plots(handles,corr_maps)
%
% % Update the correlation signal panels
% for m=1:length(handles.corr_signals)
%
% 	% Handle to the current axes
% 	axes_handle = handles.corr_signals(m);
%
% 	% Handle to the plot line
% 	line_handle = findobj(axes_handle,'Type','line','-depth',1);
%
% 	% Properties whose values we want to keep
% 	x_tick_label = get(axes_handle,'XTickLabel');
% 	y_tick_label = get(axes_handle,'YTickLabel');
% 	x_label = get(get(axes_handle,'XLabel'),'String');
% 	y_label = get(get(axes_handle,'YLabel'),'String');
% 	line_color = get(line_handle,'Color');
%
% 	% Plot the 1D correlation signal
% 	line_handle = corr_maps(m).plot_corr_signal(axes_handle);
%
% 	% Restore properties
% 	% 	set(axes_handle,'XTickLabel',x_tick_label);
% 	% 	set(axes_handle,'YTickLabel',y_tick_label);
% 	if isempty(x_tick_label), set(axes_handle,'XTickLabel',[]); end;
% 	if isempty(y_tick_label), set(axes_handle,'YTickLabel',[]); end;
% 	set(get(axes_handle,'XLabel'),'String',x_label);
% 	set(get(axes_handle,'YLabel'),'String',y_label);
% 	set(line_handle,'Color',line_color);
%
% end
%
% % Update the correlation map panels
% for m=1:length(handles.corr_maps)
%
% 	% Handle to the current axes
% 	axes_handle = handles.corr_maps(m);
%
% 	% The current position of the cursor on the correlation map. When
% 	% we call the method plot_corr_map() of the 'corr_maps' objects in
% 	% order to update the correlation map, the image object associated
% 	% with the correlation map gets deleted, along with all its children
% 	% (which includes the text object representing our cursor). So, once
% 	% the updated correlation map has been plotted, we need to create a
% 	% new cursor, place it at the same position as the old, deleted one
% 	% was, and then assign the cursor the callback function responsible
% 	% for updating its position when the correlation map is clicked on
% 	% (function update_cursor_position())
% 	cursor_position = get(handles.marker(m),'Position');
%
% 	% The current position of the colorbar
% 	colorbar_position = get(handles.colorbars(m),'Position');
%
% 	% Properties whose values we want to keep
% 	x_tick_label = get(axes_handle,'XTickLabel');
% 	y_tick_label = get(axes_handle,'YTickLabel');
% 	x_label = get(get(axes_handle,'XLabel'),'String');
% 	y_label = get(get(axes_handle,'YLabel'),'String');
%
% 	% Plot the 2D correlation map. This deletes the image object associated
% 	% with the current correlation map (along with all its children objects)
% 	% and creates a new one. We then need to update the 'handles' structure
% 	% with the new handles returned by the plot_corr_map() method
% 	[handles.corr_map_images(m),handles.colorbars(m)] = corr_maps(m).plot_corr_map(axes_handle,colorbar_position);
%
% 	% Restore properties
% 	% 	set(axes_handle,'XTickLabel',x_tick_label);
% 	% 	set(axes_handle,'YTickLabel',y_tick_label);
% 	if isempty(x_tick_label), set(axes_handle,'XTickLabel',[]); end;
% 	if isempty(y_tick_label), set(axes_handle,'YTickLabel',[]); end;
% 	set(get(axes_handle,'XLabel'),'String',x_label);
% 	set(get(axes_handle,'YLabel'),'String',y_label);
%
% 	% Ensure that the cursor is always inside the limits of the correlation
% 	% map. When we lower the offset range of the correlation map (by changing
% 	% the value of 'delta'), the cursor may fall outside of the new, smaller
% 	% correlation map. We then change the cursor position to the border of
% 	% the new correlation map. Question: do we really want this? Perhaps we
% 	% could just leave the cursor alone, as the next click on the correlation
% 	% map will update the cursor position anyways.
% 	x_lim = sort(get(axes_handle,'XLim'));
% 	y_lim = sort(get(axes_handle,'YLim'));
% 	cursor_position(1) = min(cursor_position(1),x_lim(2));
% 	cursor_position(1) = max(cursor_position(1),x_lim(1));
% 	cursor_position(2) = min(cursor_position(2),y_lim(2));
% 	cursor_position(2) = max(cursor_position(2),y_lim(1));
%
% 	% Create a new cursor and place it at the same position as the old one
% 	handles.marker(m) = text(cursor_position(1),cursor_position(2),'+','Parent',axes_handle,'Visible','On','HitTest','Off','HorizontalAlignment','Center','FontSize',18);
%
% 	% Assign a callback function to the correlation map image
% 	set(handles.corr_map_images(m),'ButtonDownFcn',{@update_cursor_position,handles.figure});
%
% end

% --------------------------------------------------------------------------------------------- %

function handles = update_corr_plots(handles,corr_maps)

% Update the correlation signal panels
for m=1:length(handles.corr_signals)
	
	% Handle to the current axes
	axes_handle = handles.corr_signals(m);
	
	% Handle to the plot line
	line_handle = findobj(axes_handle,'Type','line','-depth',1);
	
	% Update the 1D correlation map signal
	set(line_handle,'Ydata',corr_maps(m).corr_signal);
	
end

% Update the correlation map panels
for m=1:length(handles.corr_maps)
	
	% Handle to the current axes
	axes_handle = handles.corr_maps(m);
	
	% Handle to the correlation map image
	image_handle = findobj(axes_handle,'Type','image','-depth',1);
	
	% Update the correlation map's image
	set(image_handle,'CData',corr_maps(m).corr_map);
	
	% The new y-limits for both the correlation maps' axes and image
	aux = [corr_maps(m).time_offsets(1) corr_maps(m).time_offsets(end)];
	y_lim_axes  = sort(aux,2,'ascend');
	y_lim_image = sort(aux,2,'descend');
	
	% Fine tune the correlation maps' axes
	n_rows = size(corr_maps(m).corr_map,1);
	delta_y = diff(y_lim_axes)/(n_rows-1);
	y_lim_axes = y_lim_axes+[-delta_y delta_y]/2;
	
	% Set the new y-limits for the correlation map's axes and image
	set(axes_handle,'YLim',y_lim_axes);
	set(image_handle,'YData',y_lim_image);
	
	% The current position of the cursor on the correlation map
	cursor_position = get(handles.marker(m),'Position');
	
	% Check if the cursor is inside the limits of the correlation map.
	% When we lower the offset range of the correlation map (by changing
	% the value of 'delta'), the cursor may fall outside of the new axes
	% limits. If that happens, we make the cursor invisible, by setting
	% its 'Visible' property to 'off'.
	x_lim = sort(get(axes_handle,'XLim'));
	y_lim = sort(get(axes_handle,'YLim'));
	aux = (cursor_position([1 1 2 2]) > [x_lim y_lim]);
	if isequal(aux,[1 0 1 0]), status = 'On'; else status = 'Off'; end;
	set(handles.marker,'Visible',status);
	
end

% --------------------------------------------------------------------------------------------- %

function [x,x_str] = format_string(x)

% The number of decimal places used when defining the precision of 'x'
N = 3;

% Round 'x' so that it has a precision of 'N' decimal places and
% format the corresponding string accordingly
x = round(x*10^N)/10^N;
x_str = sprintf(sprintf('%%0.%gf',N),x);

% --------------------------------------------------------------------------------------------- %

function update_cursor_position(h_object,event_data,figure_handle)

% Load the shared data structure
dados = getappdata(figure_handle,'dados');

% The handles structure
handles = dados.handles;

% The handle to the correlation map where the event (mouse click) took place
corr_map_image_handle = h_object;
corr_map_handle = get(corr_map_image_handle,'Parent');

% Get the position of the clicked point
% aux = get(handles.axes_3,'CurrentPoint');
aux = get(corr_map_handle,'CurrentPoint');
mouse_position(1) = aux(1,1);
mouse_position(2) = aux(1,2);

% Get the sampling frequency from the structure 'dados'
fs = dados.corr_maps(1).rate;

% The horizontal position in signals 1 and 2 corresponding to the clicked
% position on the correlation map. Update: We used to plot the correlation
% map in a non-symmetrical way, where the second signal was the reference
% signal. Thus, the position (t,d) on the correlation map would show the
% correlation between signal_1(t-d) and signal_2(t). The correlation map
% has since been changed to be plotted in a symmetrical way and therefore
% the position (t,d) on the correlation map now shows the correlation
% between signal_1(t-d/2) and signal_2(t+d/2).
% position_1 = mouse_position(1) - mouse_position(2);
% position_2 = mouse_position(1);
position_1 = mouse_position(1) - mouse_position(2)/2;
position_2 = mouse_position(1) + mouse_position(2)/2;

% The correlation matrix
% corr_matrix = get(handles.corr_map,'CData');
corr_matrix = get(corr_map_image_handle,'CData');

% The x and y limits
% x_limits = sort(get(handles.corr_map,'XData'));
% y_limits = sort(get(handles.corr_map,'YData'));
x_limits = sort(get(corr_map_image_handle,'XData'));
y_limits = sort(get(corr_map_image_handle,'YData'));

% The mouse position and the axes limits, in number of samples. We
% need this in order to compute the index in the correlation matrix
% which corresponds to the clicked position on the correlation map
% mouse_position_samples = fix(mouse_position*fs);
mouse_position_samples = round(mouse_position*fs);
x_limits_samples = round(x_limits*fs);
y_limits_samples = round(y_limits*fs);

% When clicking on the edges of the correlation map, the position returned
% by Matlab can be slightly off the limits of the data. When converting the
% mouse position into number of samples, there may be a one-sample error,
% due to rounding. Thus, we have to check if the clicked position, in
% number of samples, is inside the data limits.
mouse_position_samples(1) = min(mouse_position_samples(1),max(x_limits_samples));
mouse_position_samples(1) = max(mouse_position_samples(1),min(x_limits_samples));
mouse_position_samples(2) = min(mouse_position_samples(2),max(y_limits_samples));
mouse_position_samples(2) = max(mouse_position_samples(2),min(y_limits_samples));

% The indice in the correlation matrix corresponding to the clicked point
indice(1) = mouse_position_samples(1);
indice(2) = y_limits_samples(2)-mouse_position_samples(2)+1;

% The value of the correlation at the clicked point
sigma = corr_matrix(indice(2),indice(1));

% % Update the text with the value of  the correlation coefficient
% texto_1 = ['corr(' num2str(position_1,4) ',' num2str(position_2,4) ') ='];
% texto_2 = num2str(sigma,4);
% set(handles.text_1,'String',texto_1);
% set(handles.text_2,'String',texto_2);

% Update the position of the marker in the correlation map
% set(handles.marker,'Position',[mouse_position 0],'Visible','On');
for m=1:length(handles.marker)
	set(handles.marker(m),'Position',[mouse_position 0],'Visible','On');
end

% Update the position of the vertical line in axes 1
line_x = [position_1 position_1];
% line_y = get(handles.axes_1,'YLim');
line_y = get(handles.signal_1,'YLim');
set(handles.line_1,'XData',line_x,'YData',line_y);

% Update the position of the vertical line in axes 2
line_x = [position_2 position_2];
% line_y = get(handles.axes_2,'YLim');
line_y = get(handles.signal_2,'YLim');
set(handles.line_2,'XData',line_x,'YData',line_y);

% --------------------------------------------------------------------------------------------- %
