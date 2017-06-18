function [handles,corr_map] = sample()

% Demonstrate the use of the instantaneous correlation algorithm and of the
% correlation map gui. The instantaneous correlation is computed between the
% tongue root of a speaker who is repeatedly saying /cop/ and the tongue tip
% of another speaker who is repeatedly saying /top/.
%
% Author: Adriano Vilela Barbosa (adriano.vilela@gmail.com)
% Copyright 2012 Adriano Vilela Barbosa


% The data file name
file_name = 'data.mat';

% Load the data file
load(file_name);

% Create TimeSignal objects from the face and forehead motion signals
signal_1 = TimeSignal(signals(1).signal,signals(1).rate);
signal_2 = TimeSignal(signals(2).signal,signals(2).rate);

% The signal names
% signal_1.name = signals(1).name;
% signal_2.name = signals(2).name;
signal_1.name = sprintf('Speaker 1\ntongue root');
signal_2.name = sprintf('Speaker 2\ntongue tip');

% ---------- The instantaneous correlation algorithm ---------- %

% The correlation map parameters
filter_type = 'ds';
eta = 0.05;
delta = round(0.5*signals(1).rate);
plot_1d = true;
interactive = true;
% plot_inputs_together = true;
plot_inputs_together = false;

% Plot the correlation map
[handles,corr_map] = correlation_map_gui(signal_1,signal_2,[],filter_type,eta,delta,plot_1d,interactive,plot_inputs_together);

% ---------------------------------------------------------------------------- %
