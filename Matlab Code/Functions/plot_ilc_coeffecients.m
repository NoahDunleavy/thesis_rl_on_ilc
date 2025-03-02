function [alpha_error_fig, beta_error_fig, alpha_prog_fig, beta_prog_fig] = plot_ilc_coeffecients(figure_name, graph_title, ILC_Trial, to_plot, beta_star, betas_to_plot, alphas_to_plot, save_path)
%Note: beta star error is sometimes illogical
%Plot out the progression of coeffecients in an ILC problem described in a
%basis space
%Inputs:
    %figure_name: string - name of the figure that is opened, or precursor to saved image
    %graph_title: string - displayed title header on plot
    %ILC_Trial: structure - containing iterations of coeffecients
    %to_plot: scalar/vector - indicate which trials or how many to plot
    %beta_star: vector - sometimes illogical, but when there is a defined beta weights, use to compute the error
    %betas_to_plot: scalar/vector - indicate which beta values (or how many) to plot
    %alphas_to_plot: scalar/vector - see above
    %save_path: string - where to save the generated plots to
%Outputs:
    %alpha_error_fig: figure - handle to alpha errors
    %beta_error_fig: figure - handle to beta erros
    %alpha_prog_fig: figure - handle to showing alpha evolve
    %beta_prog_fig: figure - handle to showing betas evolve

max_coef = 5;   %default maximum # of coeffeceints to plot
if nargin < 8   %no save path
    save_path = -1;
end
if nargin < 7 %if no alpha to plot
    alphas_to_plot = max_coef;
end
if nargin < 6
    betas_to_plot = max_coef;
end
if nargin < 5
    beta_star = -1;    %ensure we have a mtrix for mathing
end
beta_star_passed = true;
if isscalar(beta_star)
    if(beta_star == -1)
        beta_star_passed= false;
    end
end

num_trials = length(ILC_Trial);
marker_color_scale = 0.75;  %how to scale the color difference from line to point
marker_size = 20;
subtitle_size = getappdata(groot, 'DefaultSubtitleFontSize');

if isscalar(to_plot)    %check if we are passing in a count or specific trials to plot
    trials_to_plot = floor(linspace(1, num_trials, to_plot));
    trials_to_plot(end) = num_trials;
else
    trials_to_plot = to_plot;
end

%Process ILC Strcture
num_betas = length(ILC_Trial(1).betas);
num_alphas = length(ILC_Trial(1).alphas);
num_plot = length(trials_to_plot);

%Determine which coeffecients to plot
if isscalar(alphas_to_plot)
    if(alphas_to_plot == -1)    %default plot
        alphas_to_plot = max_coef;
    end
    alphas_to_plot = floor(linspace(1, num_alphas, alphas_to_plot)); %space out the plots 
    alphas_to_plot(end) = num_alphas; %ensure the last one is shown
end
if isscalar(betas_to_plot)
    if(betas_to_plot == -1)    %default plot
        betas_to_plot = max_coef;
    end
    betas_to_plot = floor(linspace(1, num_alphas, betas_to_plot)); %space out the plots 
    betas_to_plot(end) = num_betas; %ensure the last one is shown
end

%Progression of errors on coeffecients
alpha_error_norm = zeros(num_trials, 1);
beta_error_norm = zeros(num_trials, 1);
for ndx = 1:num_trials
    alpha_error_norm(ndx) = norm(ILC_Trial(ndx).alpha_error) / num_alphas;
    if beta_star_passed
        beta_error_norm(ndx) = norm(beta_star - ILC_Trial(ndx).betas) / num_betas;
    end
end

alpha_error_fig = figure('Name', sprintf('%s - Alpha Error Norm Progression', figure_name));
plot(1:num_trials, alpha_error_norm, 'r*', 'LineStyle',':');
title('Alpha Error', graph_title);
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Trial Number')
ylabel('Normalized Amplitude')

if beta_star_passed %if we were given a beta_star
    beta_error_fig = figure('Name', sprintf('%s - Beta Error Norm Progression', figure_name));
    plot(1:num_trials, beta_error_norm, 'r*', 'LineStyle',':');
    title('Beta Error', graph_title)
    subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
    xlabel('Trial Number')
    ylabel('Normalized Amplitude')
end

%Transfer structure to array for plotting
%For progression smoothing
full_betas = zeros(num_betas, num_trials);
full_alphas = zeros(num_alphas, num_trials);
for ndx = 1:num_trials
    full_betas(:, ndx) = ILC_Trial(ndx).betas;
    full_alphas(:, ndx) = ILC_Trial(ndx).alphas;
end

%Progression of Alphas
alpha_prog_fig = figure('Name', sprintf('%s - Alpha Coeffecient Progression', figure_name));
legend_list = "[REMOVE ME]";
hold on;
for alpha = alphas_to_plot
    temp_plot = plot(0:(num_trials-1), full_alphas(alpha, :), 'MarkerSize', marker_size, 'LineStyle', '-', 'Marker', 'o', 'MarkerIndices', trials_to_plot);
    temp_plot.MarkerFaceColor = temp_plot.Color * marker_color_scale;%make the marker colors slightly darker than lines
    temp_plot.MarkerEdgeColor = temp_plot.Color * marker_color_scale;
    legend_list(end + 1) = sprintf('Alpha %.d', alpha);
end
hold off;

title('Alpha Coeffecients')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Trial Number')
ylabel('Alpha Value')
xticks(trials_to_plot - 1);  
legend(legend_list(2:end))

%Progression of Beats
beta_prog_fig = figure('Name', sprintf('%s - Beta Coeffecient Progression', figure_name));
legend_list = "[REMOVE ME]";
hold on;
for beta = betas_to_plot
    temp_plot = plot(0:(num_trials-1), full_betas(beta, :), 'MarkerSize', marker_size, 'LineStyle', '-', 'Marker', 'o', 'MarkerIndices', trials_to_plot);
    temp_plot.MarkerFaceColor = temp_plot.Color * marker_color_scale; %make the marker colors slightly darker than lines
    temp_plot.MarkerEdgeColor = temp_plot.Color * marker_color_scale;
    legend_list(end + 1) = sprintf('Beta %.d', beta);
end
hold off;

title('Beta Coeffecients')
subtitle(sprintf('%s', graph_title), 'FontSize', subtitle_size)
xlabel('Trial Number')
ylabel('Beta Value')
xticks(trials_to_plot - 1);  
legend(legend_list(2:end))

%Save images if possible
if exist('save_path', 'var')    %if a save path was provided
    save_figure(save_path, [alpha_error_fig, alpha_prog_fig, beta_prog_fig]); %save all the figures to the path
    if beta_star_passed
        save_figure(save_path, beta_error_fig);
    end
end

end
