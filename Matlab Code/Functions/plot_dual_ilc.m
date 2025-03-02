function [error_fig, out1_fig, out2_fig, dual_out_fig, input1_fig, input2_fig] = plot_dual_ilc(figure_name, graph_title, ILC_Trial, goal_output, goal_input, to_plot, save_path)
%Plot the ILC Structure information, showing the input and output history
%for a two-spring-mass system through ILC trials
%Inputs:
    %figure_name: string - name of the figure that is opened, or precursor to saved image
    %graph_title: string - displayed title header on plot
    %ILC_Trial: structure - containing iterations of coeffecients
    %goal_output: vector - y*, what we want to generate
    %goal_input: vector - u*, what input gets us there
    %to_plot: scalar/vector - indicate which trials or how many to plot
    %save_path: string - where to save the generated plots to
%Outputs:
    %error_fig: figure - handle to the error progression
    %out1_fig: figure - handle to position of mass 1
    %out2_fig: figure - handle to position of mass 2
    %dual_out_fig: figure - handle to figure showing a 'shaped' output
    %input1_fig: figure - handle to input 1
    %input2_fig: figure - handle to input 2

%legend_alpha = 0.5;
subtitle_size = getappdata(groot, 'DefaultSubtitleFontSize');

if nargin < 7   %if no file save path
    save_path = -1;
end
if nargin < 6
    to_plot = 5;    %default to plotting 5 trials
end
if nargin < 5
    goal_input = -1;
end
if nargin <4
    goal_output = -1;
end


%Process ILC Strcture
num_trials = length(ILC_Trial);
if isscalar(to_plot)    %check if we are passing in a count or specific trials to plot
    trials_to_plot = floor(linspace(1, num_trials, to_plot));
    trials_to_plot(end) = num_trials;
else
    trials_to_plot = to_plot;
end

%There is potential to clean up the code here, but readability is more
%important
num_inputs = length(ILC_Trial(1).input);
num_outputs = length(ILC_Trial(1).output);

out1_ndx = 1:2:num_outputs;
out2_ndx = 2:2:num_outputs;

in1_ndx = 1:2:num_inputs;
in2_ndx = 2:2:num_inputs;

%Error Progression
error_progression = zeros(num_trials, 1);
for ndx = 1:num_trials       %convert structure to array
    error_progression(ndx) = norm(ILC_Trial(ndx).output_error)/num_outputs; %normalize the errors
end
error_fig = figure('Name', sprintf('%s - Error Progression', figure_name));
plot(1:num_trials, error_progression, 'r*', 'LineStyle',':');
title('Error Magnitude');
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Trial Number')
ylabel('Normalized Amplitude')

%Mass 1 Position
out1_fig = figure('Name', sprintf('%s - Mass 1 Position', figure_name));
legend_list = "[REMOVE ME]";    %teach it we're building strings
hold on;
for ndx = 1:length(trials_to_plot)
    trial_num = trials_to_plot(ndx);    %which trial of the process is being plotted
    stairs(0:(num_outputs/2-1), ILC_Trial(trial_num).output(out1_ndx));
    legend_list(end + 1) = sprintf('Trial %.d', trial_num);
end
if ~(isscalar(goal_output) && (goal_output == -1))
    stairs(1:((num_outputs-1)/2), goal_output(out1_ndx(1:(end-1))), 'Color', [1, 0, 0], 'LineStyle', '--');
    legend_list(end + 1) = 'Goal Output';
end
hold off;
legend(legend_list(2:end));%, "BackgroundAlpha", legend_alpha);
title('Mass 1 Position')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Position (m)')

%Mass 2 Position
out2_fig = figure('Name', sprintf('%s - Mass 2 Position', figure_name));
legend_list = "[REMOVE ME]";    %reset legend list
hold on;
for ndx = 1:length(trials_to_plot)
    trial_num = trials_to_plot(ndx);    %which trial of the process is being plotted
    stairs(0:(num_outputs/2-1), ILC_Trial(trial_num).output(out2_ndx));
    legend_list(end + 1) = sprintf('Trial %.d', trial_num);
end

if ~(isscalar(goal_output) && (goal_output == -1))
    stairs(1:((num_outputs-1)/2), goal_output(out2_ndx(1:(end-1))), 'Color', [1, 0, 0], 'LineStyle', '--');
    legend_list(end + 1) = 'Goal Output';
end
hold off;
legend(legend_list(2:end));%, "BackgroundAlpha", legend_alpha);
title('Mass 2 Position')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Position (m)')

%Combined Outputs (for fun 2D shapes)
dual_out_fig = figure('Name', sprintf('%s - Shaped Output', figure_name));
legend_list = "[REMOVE ME]";    %reset legend list
hold on;
for ndx = 1:length(trials_to_plot)
    trial_num = trials_to_plot(ndx);    %which trial of the process is being plotted
    plot(ILC_Trial(trial_num).output(out1_ndx), ILC_Trial(trial_num).output(out2_ndx));
    legend_list(end + 1) = sprintf('Trial %.d', trial_num);
end
if ~(isscalar(goal_output) && (goal_output == -1))
    plot(goal_output(out1_ndx(1:(end-1))), goal_output(out2_ndx(1:(end-1))), 'Color', [1, 0, 0], 'LineStyle', '--');
    legend_list(end + 1) = 'Goal Output';
    %Add a start/stop position indicator
    scatter(goal_output(1), goal_output(2), 'filled', '>', 'MarkerFaceColor', [0, 1, 0])
    scatter(goal_output(end - 1), goal_output(end), 'filled', 'hexagram', 'MarkerFaceColor', [1, 0, 0])
    legend([legend_list(2:end), 'Start', 'Stop']);%, "BackgroundAlpha", legend_alpha);
end
hold off;
title('Shaped Outputs')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Mass 1 Position (m)')    
ylabel('Mass 2 Position (m)')
axis equal %shapes should be equally scaled


%Input 1
input1_fig = figure('Name', sprintf('%s - Input 1', figure_name));
legend_list = "[REMOVE ME]";    %reset legend list
hold on;
for ndx = 1:length(trials_to_plot)
    trial_num = trials_to_plot(ndx);    %which trial of the process is being plotted
    stairs(0:(num_inputs/2-1), ILC_Trial(trial_num).input(in1_ndx));
    legend_list(end + 1) = sprintf('Trial %.d', trial_num);
end
if ~(isscalar(goal_input) && (goal_input == -1))
    stairs(0:((num_inputs/2)-1), goal_input(in1_ndx), 'Color', [1, 0, 0], 'LineStyle', '--');
    legend_list(end + 1) = 'Goal Input';
end
hold off;
legend(legend_list(2:end));%, "BackgroundAlpha", legend_alpha);
title('Input 1')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Input (N)')

%Input 2
input2_fig = figure('Name', sprintf('%s - Input 2', figure_name));
legend_list = "[REMOVE ME]";    %reset legend list
hold on;
for ndx = 1:length(trials_to_plot)
    trial_num = trials_to_plot(ndx);    %which trial of the process is being plotted
    stairs(0:(num_inputs/2-1), ILC_Trial(trial_num).input(in2_ndx));
    legend_list(end + 1) = sprintf('Trial %.d', trial_num);
end
if ~(isscalar(goal_input) && (goal_input == -1))
    stairs(0:((num_inputs/2)-1), goal_input(in2_ndx), 'Color', [1, 0, 0], 'LineStyle', '--');
    legend_list(end + 1) = 'Goal Input';
end
hold off;
legend(legend_list(2:end));%, "BackgroundAlpha", legend_alpha);
title('Input 2')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Input (N)')

%Save images if possible
if exist('save_path', 'var')    %if a save path was provided
    save_figure(save_path, [error_fig, out1_fig, out2_fig, input1_fig, input2_fig]); %save all the figures to the path
    save_figure(save_path, dual_out_fig, -1, true); %this one is shapede
end

end