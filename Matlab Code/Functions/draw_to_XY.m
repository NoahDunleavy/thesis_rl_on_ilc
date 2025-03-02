function [x, y] = draw_to_XY(resolution)
    %Written by ChatGPT with prompt: 'I would like a matlab function that 
    %opens up a figure that I can draw on with my mouse, and then it returns 
    %two arrays of the x positions and y positiions of the thing I drew. I 
    %would like to set the resoltuon (aka number of xy pairs) in the 
    %function call'

    % draw_to_XY - Opens a figure to draw with the mouse and returns x, y coordinates.
    %
    % Syntax: [x, y] = draw_with_mouse(resolution)
    %
    % Inputs:
    %   resolution - Number of points to return in the output arrays.
    %
    % Outputs:
    %   x - Array of x positions of the drawn curve.
    %   y - Array of y positions of the drawn curve.

    if nargin < 1
        resolution = 100; % Default resolution if not specified
    end

    % Create a new figure
    figure;
    hold on;
    grid on
    axis([0 1 0 1]); % Set axis limits (can be adjusted as needed)
    title('Draw with your mouse (click and drag). Press Enter when done.');
    set(gca, 'Position', [0.1, 0.1, 0.8, 0.8]); % Adjust axis position for aesthetics
    set(gcf, 'WindowButtonMotionFcn', @mouse_move_callback); % Enable mouse motion tracking
    set(gcf, 'WindowButtonDownFcn', @mouse_down_callback); % Start drawing on mouse press
    set(gcf, 'WindowButtonUpFcn', @mouse_up_callback); % Stop drawing on mouse release
    pan off; zoom off; % Disable default interactions
    
    % Variables to store points
    drawing = false; % Indicates whether the mouse is pressed
    points_x = [];
    points_y = [];

    % Callback to track mouse movements while pressed
    function mouse_move_callback(~, ~)
        if drawing
            current_point = get(gca, 'CurrentPoint');
            points_x(end + 1) = current_point(1, 1);
            points_y(end + 1) = current_point(1, 2);
            plot(points_x, points_y, 'b-');
            drawnow;
        end
    end

    % Callback to start drawing
    function mouse_down_callback(~, ~)
        drawing = true;
    end

    % Callback to stop drawing
    function mouse_up_callback(~, ~)
        drawing = false;
    end

    % Wait for the user to press Enter
    pause;
    close(gcf); % Close the figure

    % Interpolate to match the desired resolution
    t_raw = linspace(0, 1, length(points_x));
    t_interp = linspace(0, 1, resolution);

    x = interp1(t_raw, points_x, t_interp, 'linear');
    y = interp1(t_raw, points_y, t_interp, 'linear');
end
