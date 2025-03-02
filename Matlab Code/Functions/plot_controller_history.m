function [figure_list] = plot_controller_history(figure_name, graph_title, controllers, inputs_of_interest, states_of_interest, save_path, decoupled)
%Plot the progression of a given controller parameters for indicated inputs
%Inputs:
    %figure_name: string - name of the figure that is opened, or precursor to saved image
    %graph_title: string - displayed title header on plot
    %controllers: matrix - r x n x l matrix, showing the learning of r controllers from n states over l iterations
    %inputs_of_interest: scalar/vector - for when there are lots of inputs, which ones to plot
    %states_of_interest: scalar/vector - see above
    %save_path: string - where to save the generated plots to
    %decoupled: bool - indicates whether or not it was learned via input decoupling (only show trials where input changes)
%Outputs:
    %figure_list: list - figure handles generated

%Plot settings
marker_color_scale = 0.75;  %how to scale the color difference from line to point
marker_size = 20;
subtitle_size = getappdata(groot, 'DefaultSubtitleFontSize'); %default set in fig properties

num_inputs = height(controllers);
num_states = width(controllers);

%Assign defaultss
if nargin < 7   %if flag of whether or not it was input decoupled
    decoupled = false;
end
if nargin < 6   
    save_path = -1;
end
if nargin < 5   %if no states index passed
    states_of_interest = floor(linspace(1, num_states, min(num_states, 4)));    %do not overclutter sith state points
end
if nargin < 4   %if no passed input indexes
    inputs_of_interest = [1, num_inputs];  %default to showing 2 of them (sometimes there are a lot)
end

%Shortcuts
if (states_of_interest == -1)   %if -1, plot all components
    states_of_interest = 1:num_states;
elseif (isscalar(states_of_interest))
    states_of_interest = floor(linspace(1, num_states, states_of_interest));
    states_of_interest(end) = num_states;
end

if (inputs_of_interest == -1)
    inputs_of_interest = 1:num_inputs;
elseif (isscalar(inputs_of_interest))
    inputs_of_interest = floor(linspace(1, num_inputs, inputs_of_interest));
    inputs_of_interest(end) = num_inputs;
end

figure_list = [];
for input_num = inputs_of_interest
    figure_list(end + 1) = figure('Name', sprintf('%s - Input %.d Controller Weights', figure_name, input_num));
    
    if decoupled
        trial_ndx = [1, (input_num+1):num_inputs:size(controllers, 3)];    %if input decoupled, alternate to only show points when updated
    else
        trial_ndx = 1:size(controllers, 3);
    end

    legend_list = "[REMOVE ME]";
    hold on;
    for state = states_of_interest
        temp_plot = plot((trial_ndx-1), squeeze(controllers(input_num, state, trial_ndx)), 'MarkerSize', marker_size, 'Marker', 'o');
        temp_plot.MarkerFaceColor = temp_plot.Color * marker_color_scale;%make the marker colors slightly darker than lines
        temp_plot.MarkerEdgeColor = temp_plot.Color * marker_color_scale;
        legend_list(end + 1) = sprintf('Weight on State %.d', state);
    end
    hold off;

    legend(legend_list(2:end), "BackgroundAlpha", 0.5)
    title(sprintf('Input %.d Controller Weights', input_num));
    subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
    xlabel('Controller Number')
    ylabel('Controller Weight')
    xticks(trial_ndx-1);   


end

%Save images if possible
if exist('save_path', 'var')    %if a save path was provided
    for fig_num = figure_list
        save_figure(save_path, figure(fig_num)); %save all the figures to the path
    end
end


end