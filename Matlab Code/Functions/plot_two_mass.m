function [mass1_fig, mass2_fig, input1_fig, input2_fig] = plot_two_mass(figure_name, graph_title, outputs, inputs, save_path)
%For the dual spring-mass system repeatedly used, plot positions and inputs
%Inputs:
    %figure_name: string - name of the figure that is opened, or precursor to saved image
    %graph_title: string - displayed title header on plot
    %outputs: matrix - 2 x k matrix, where row 1 is mass 1 position, row 2 is mass 2 position
    %inputs: matrix - 2 x k, row 1 is input1 and row 2 is input 2
    %save_path: string - where to save the generated plots to
%Outputs:
    %mass1_fig: figure - handle to position of mass 1
    %mass2_fig: figure - handle to position of mass 2
    %input1_fig: figure - handle to input 1
    %input2_fig: figure - handle to input 2


if (height(outputs) > width(outputs))    %Ensure wide
    outputs = outputs';
end

if (height(inputs) > width(inputs))    %Ensure wide
    inputs = inputs';
end

num_samples = width(outputs);
subtitle_size = getappdata(groot, 'DefaultSubtitleFontSize');

%Output Style
output_color = [0    0.4471    0.7412];
output_size = 1.2;
output_style = '-';

%Mass 1 Position
mass1_fig = figure('Name', sprintf('%s - Mass 1 Position', figure_name));
%stairs(0:(num_samples-1), outputs(1, :), 'Color', output_color, 'LineStyle', output_style, 'LineWidth', output_size);  
stairs(0:(num_samples-1), outputs(1, :), 'Color', output_color);
title('Mass 1 Position');
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Position (m)')

%Mass 2 Position
mass2_fig = figure('Name', sprintf('%s - Mass 2 Position', figure_name));
stairs(0:(num_samples-1), outputs(2, :), 'Color', output_color);  
title('Mass 2 Position')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Position (m)')

%Input Style
input_color = [1, 0, 1];
input_size = 1.2;
input_style = '-';

%Input 1
input1_fig = figure('Name', sprintf('%s - Input 1 Magnitude', figure_name));
stairs(0:(num_samples-1), inputs(1, :), 'Color', input_color);  
title('Input 1')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Force (N)')

%Input 2
input2_fig = figure('Name', sprintf('%s - Input 2 Magnitude', figure_name));
stairs(0:(num_samples-1), inputs(2, :), 'Color', input_color);  
title('Input 2')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Sample Number (k)')    
ylabel('Force (N)')

%Save images if possible
if exist('save_path', 'var')    %if a save path was provided
    save_figure(save_path, [mass1_fig, mass2_fig, input1_fig, input2_fig]); %save all the figures to the path
end

end