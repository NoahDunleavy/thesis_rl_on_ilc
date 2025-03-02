function [fig_pole] = plot_pole_placement(figure_name, graph_title, pole_locations, goal_poles, save_path)
%Render the pole placements (or more accurately, plot eigen values)
%Inputs:
    %figure_name: string - name of the figure that is opened, or precursor to saved image
    %graph_title: string - displayed title header on plot
    %pole_locations: vector - location of poles
    %goal_poles: vector - where we hoped/wanted poles (not rendered if excluded)
    %save_path: string - location to save the figure
%Outputs:
    %fig_pole: figure - handle to showing poles

if ((nargin < 4) || (all(goal_poles == -1)))  %if there is no goal poles
    goal_poles = false;
end

fig_pole = figure('Name', figure_name);
circle_resolution = 200;    %how many points to make up the circle
circle_domain = linspace(0, 2*pi, circle_resolution);    %trace out a full period
unit_x = cos(circle_domain);
unit_y = sin(circle_domain);
plot(unit_x, unit_y);
hold on;
scatter(real(pole_locations), imag(pole_locations), 'filled', 'o', 'LineWidth', 3, 'SizeData', 200);

xline(0)
yline(0)
if (goal_poles == false) %yes I know usually you can just do ~goal_poles, but if this is a matrix that breaks it
    legend('Unit Circle', 'Placed Poles')
else
    scatter(real(goal_poles), imag(goal_poles), 'o', 'SizeData', 200);
    legend('Unit Circle', 'Pole Locations', '', '', 'Goal Poles')   %skip the x and y line
end
hold off;

xlabel('Real')
ylabel('Imaginary')
title('Pole Locations')
subtitle(graph_title, 'FontSize', getappdata(groot, 'DefaultSubtitleFontSize'));
axis equal  %ensure the circle looks like a circle


%Save images if possible
if exist('save_path', 'var')    %if a save path was provided
    save_figure(save_path, fig_pole, -1, true); %save all the figures to the path, marking it as shaped
end

end