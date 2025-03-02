function [file_path] = save_figure(path, fig_list, name_overwrite, shaped)
%Save figures to specified path following naming convention
%Inputs:
    %path: string - where to save the file to, 
    %fig_list: list - figure handles
    %name_overwrite: string - if you wish to overwrite the naming convention
    %shaped: bool - whether or not it is shaped (should be square)
%Outputs:
    %file_path: string - file save location

%Set save settings
title_size = 45;
subtitle_size = 30;
axes_size = 40;
tick_size = 35;
legend_size = 30;

%Default condition
if ~exist('name_overwrite', 'var')
    name_overwrite = -1;
end
if ~exist('shaped', 'var')
    shaped = false;
end

if path == -1   %condition for if we are not updating the figures
    return
end

%Dimension params
aspect_ratio = 16/12;
figure_width = 16;
figure_height = figure_width / aspect_ratio;
if shaped
    figure_width = figure_height;
end

for fig = fig_list
    %Determine what name to give the file
    if name_overwrite == -1
        file_name = fig.Name;
    else
        file_name = name_overwrite;
    end

    %Set figure sizes
    fig.Units = 'inches';
    fig.Position = [0, 0, figure_width, figure_height];

    % Ensure the axes is square
    ax = findobj(fig, 'Type', 'axes');
    for a = ax'
        %Set font sizes
        a.FontSize = tick_size; %tick size - do this one first because otherwise it overwrites others
        title(a, a.Title.String, 'FontSize', title_size, 'FontWeight', 'bold');
        xlabel(a, a.XLabel.String, 'FontSize', axes_size);
        ylabel(a, a.YLabel.String, 'FontSize', axes_size);

        %adjust subtitles
        if isprop(a, 'Subtitle') && ~isempty(a.Subtitle.String)
            subtitle(a, a.Subtitle.String, 'FontSize', subtitle_size, 'FontWeight', 'normal');
        end

        %adjust legend
        lgd = findobj(fig, 'Type', 'legend');
        for l = lgd'
            l.FontSize = legend_size;
            l.Location = 'northeast';
        end

        set(a, 'PlotBoxAspectRatio', [figure_width, figure_height, 1]); %ensure plot matched to our desired scales
        if shaped
            set(a, 'DataAspectRatio', [1, 1, 1]); %force square
            axis(a, 'equal'); %dont distort scaling
        else
            set(a, 'DataAspectRatioMode', 'auto');
            axis(a, 'tight'); % Ensure scaling is not distorted
        end
        %Fix offset caused by labels
        insets = get(a, 'TightInset'); %space needed for labels
        left_offset = insets(1); %padding due to ylabel
        bottom_offset = insets(2); %padding due to xlabel
        right_offset = insets(3); %padding due to xlabel
        top_offset = insets(4); %padding due to title

        % Adjust outer position to fit tightly around the plot content
        % Ensure no unnecessary margins are included
        a.OuterPosition = [left_offset/figure_width, bottom_offset/figure_height, (figure_width - left_offset - right_offset)/figure_width, (figure_height - top_offset - bottom_offset)/figure_height];
    end

    %setup space such that image makes room for outer elements (title,
    %axes, etc)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0, 0, figure_width, figure_height]; %make saved figure match rendered image
    fig.PaperSize = [figure_width, figure_height]; %make paper and fig match
    drawnow;

    file_path = sprintf('%s\\%s.pdf', path, file_name);
    print(fig, file_path, '-dpdf', '-r300');  % pdf for LaTeX

end

end