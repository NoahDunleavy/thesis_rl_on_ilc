%A General Intro to State Space and Modern Control Theory
%Noah Dunleavy
%Honors Thesis under direction of Professor Minh Q. Phan
%Thayer School of Engineering
clc; clear;
addpath('Saved Data', 'Functions');
setDefaultFigProp();
plot_all = false;    %most run time goes into generating figs
update_file_path = -1; %'C:\Users\noahd\OneDrive\Desktop\Thesis\Thesis Images\General Intro';  %set this to the save location if want to update figures, or -1 if not
if update_file_path ~= -1
    keyboard    %ensure we have to be very concious about ever updating all the images
end

%% Description
%{
This intro code provides the basics of 'known' status in the field. Tot hose familiar with the Thayer School of Engineering's curriculum, it is a recap of ENGS145 - Modern Control Theory.

We begin with system formulation and representation in the matrix form, and how we resolve the difference between the continuous nature of the world, and the discrete ability of computers. Here, we use the 'Zero-Order-Hold' approach. Next the idea of pole placement is introduced, and it is demonstrated how the further from the origin the poles are, the longer control takes. A deadbeat controller is used to highlight this, which is the time-optimal solution for any system.

One will note that the deadbeat controller, while time optimal, requires significant control effort that may not be realsitic or safe for a real system. That leads us to our introduction of the Linear Quadratic Regulation (LQR) controller, which minimizes a cost function composed of both system inputs and states. We show the difference in responses when input has two different relative strengths when compared to the state. 

Next we introduce the Iterative Learning Control (ILC) problem and a controller, and show that it can learn to generate any output (so long as permitted by the physical characteristics of the system). 

All of the previous examples rely on perfect knowledge not typically possible in the real world (although we can approximate and identiy with System Identification Methods), the process of Reinforcement Learning (RL) is shown via the Policy Iteration method. It can be, and is, shown to find the LQR controller as defined by its cost function. 

Finally, the idea of basis functions is introduced briefly to lay the
groundwork for subsequent explorations.
%}

%% System Creation

%Scalar Representation
m1 = 1;     %Mass of the block, [kg]
m2 = 0.5;
num_masses = 2;
k1 = 100;   %Spring constant, [N/m]
k2 = 200;
c1 = 1;  %Dampning Coeffeceint, [Ns/m]
c2 = 0.5;

%Matrix Formulations
M = diag([m1, m2]);   %mass matrix
K = [k1 + k2, -k2; -k2, k2];  %stiffness matrix
C_damp = [c1 + c2, -c2; -c2, c2] ; %Dampning matrix

%Continuous Statespace Formulations
Ac = [zeros(num_masses), eye(num_masses); -M^-1 * K, -M^-1 * C_damp];
Bc = [zeros(num_masses); M^-1];
C = [1, 0, 0, 0; 0, 1, 0, 0];   %Monitor the block positions as outputs
D = [0, 0; 0, 0];    %

%Dimensional Variables
num_states = height(Ac);    %Ac is n x n
num_inputs = width(Bc);     %Bc is n x r
num_outputs = height(C);    %C  is m x n

%Ensure against Nyquist
eigens = eig(Ac);                   %Ac captures the modes of a given system
natural_freqs = imag(eigens)/(2*pi);%Get all the natural frequencies in Hz (convert from rad/sec)
max_freq = max(natural_freqs);      %The highest frequency is the one we must cater to
nyquist_freq = 2 * max_freq;        %Minimum sampling frequency to avoid nyquist sampling
nyquist_rate = 1 / nyquist_freq;    %Sample frequency to period

dt = 0.01;              %Set the sampling rate we will be using
if (dt > nyquist_rate)  %Ensure that the true sampling period is below the nyquist specifications
    dt = nyquist_rate / 10; %Set the dt to one that is sufficient, best practice is to scale by a factor of 2-10 beyond nyquist
    sprintf('Selected sample rate was insufficient, setting dt = %.2e seconds', dt)
end

%Check Controllability
if ~is_controllable(Ac, Bc) %always check it is even possible to control
    printf('The Continuous System is not Controllable!')
end

%Convert to Discrete
[A, B] = c2d(Ac, Bc, dt);   %Performs the necessary discretization operations - verify this for yourself

%Set initial conditions
initial_displacement = [1; 0];   %displace the cooresponding blocks initially by this amount, [meters]
initial_velocity = [0; 2];       %Start cooresponding blocks with this velocity [m/s]
x0 = [initial_displacement; initial_velocity];  %The top half sets position ICs, and the bottom half sets velocity

%Simulation Duartion
max_samples = 500;           %set the maximum number of steps out to simulate (note: max k will be -1, since we start at k=0)
max_time = max_samples * dt;  %Important distinction: sample STEPS is not the same as time

save('Saved Data\thesis_system.mat', 'A', 'B', 'C', 'D', 'x0', 'num_states', 'num_outputs', 'num_inputs'); %save the system for consistency across all following code
%% Open-Loop Simulation
cont_vs_disc_resolution = 1;   %Additional resolution to give continuous sinulations for smoother plots, and so the evident differences can be seen between continious and discrete
%Setting to 1 speeds up rendering process, and can just use matlabs
%'plot' vs stairs to highlight

open_time_continuous = linspace(0, max_time, max_samples * cont_vs_disc_resolution)';%Generate the time scale at the higher resolution
continuous_system = ss(Ac, Bc, C, D);   %Turn the continuous matricies into a system for simulation
continuous_open_outputs = lsim(continuous_system, zeros(height(open_time_continuous), num_inputs), open_time_continuous, x0);
[discrete_open_outputs, discrete_open_states] = dlsim(A, B, C, D, zeros(max_samples, num_inputs), x0);     %Generate x(0) -> x(max_samples - 1), and cooresponding outputs (y)

%% Render Open Loop Responses
if plot_all
    %Block 1 Position
    fig = figure('Name', 'Open-loop Mass 1 Position');
    plot(open_time_continuous, continuous_open_outputs(:, 1));

    title('Mass 1 Position');
    subtitle('Open-Loop', 'FontSize', getappdata(groot, 'DefaultSubtitleFontSize'));
    xlabel('Time (sec)')    %only point where we will be using time as the x axis
    ylabel('Position (m)')
    save_figure(update_file_path, fig, 'Continuous Open-Loop - Mass 1'); %save just the continuous Model
    %Add in discrete element
    hold on;
    stairs(dt * (0:(max_samples-1)), discrete_open_outputs(:, 1));  %note the dt scaling such that axises are the same
    hold off;
    legend('Continuous', 'Discrete')
    save_figure(update_file_path, fig);
    
    %Mass 2 position
    fig = figure('Name', 'Open-loop Mass 2 Position');
    plot(open_time_continuous, continuous_open_outputs(:, 2));
    title('Mass 2 Position');
    subtitle('Open-Loop', 'FontSize', getappdata(groot, 'DefaultSubtitleFontSize'));
    xlabel('Time (sec)')
    ylabel('Position (m)')
    save_figure(update_file_path, fig, 'Continuous Open-Loop - Mass 2'); %save just the continuous Model
    %Add in discrete element
    hold on;
    stairs(dt * (0:(max_samples-1)), discrete_open_outputs(:, 2));
    hold off;
    legend('Continuous', 'Discrete')
    save_figure(update_file_path, fig);

    %Zoomed in to show the models match
    num_zoomed_in = 11;
    fig = figure('Name', 'Zoomed-in Open-loop Mass 1 Position');
    plot(open_time_continuous(1:num_zoomed_in), continuous_open_outputs(1:num_zoomed_in, 1));
    hold on;
    stairs(dt * (0:(num_zoomed_in-1)), discrete_open_outputs(1:num_zoomed_in, 1));  %note the dt scaling such that axises are the same
    scatter(dt * (0:(num_zoomed_in-1)), discrete_open_outputs(1:num_zoomed_in, 1), 'filled', '*', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'SizeData', 90);
    hold off;
    legend('Continuous', 'Discrete')
    title('Mass 1 Position');
    subtitle('Open-Loop', 'FontSize', getappdata(groot, 'DefaultSubtitleFontSize'));
    xlabel('Time (sec)')    %only point where we will be using time as the x axis
    ylabel('Position (m)')
    xlim([0, dt*(num_zoomed_in-1)])
    save_figure(update_file_path, fig);

    fig = figure('Name', 'Zoomed-in Open-loop Mass 2 Position');
    plot(open_time_continuous(1:num_zoomed_in), continuous_open_outputs(1:num_zoomed_in, 2));
    hold on;
    stairs(dt * (0:(num_zoomed_in-1)), discrete_open_outputs(1:num_zoomed_in, 2));  %note the dt scaling such that axises are the same
    scatter(dt * (0:(num_zoomed_in-1)), discrete_open_outputs(1:num_zoomed_in, 2), 'filled', '*', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'SizeData', 90);
    hold off;
    legend('Continuous', 'Discrete')
    title('Mass 2 Position');
    subtitle('Open-Loop', 'FontSize', getappdata(groot, 'DefaultSubtitleFontSize'));
    xlabel('Time (sec)')    %only point where we will be using time as the x axis
    ylabel('Position (m)')
    xlim([0, dt*(num_zoomed_in-1)])
    save_figure(update_file_path, fig);
end

%% Simple Pole Placement Controller
pole_locations = [0.5 + 0.5i, 0.5 - 0.5i, -0.7 + 0.1i, -0.7 - 0.1i];   %if imaginary, must be complex conjugates of eachother
placed_controller = -place(A, B, pole_locations);
verify_poles = eig(A + B * placed_controller);
if any(abs(verify_poles) >= 1)
    printf('A Pole was placed outside of the unit circle, unstable controller!')
end

%Visualize the poles
if plot_all
    plot_pole_placement('Simple Pole Placement', 'Simple Pole Placement', verify_poles, pole_locations, update_file_path);
end
closed_pole_samples = 60;
pole_states = zeros(num_states, closed_pole_samples);
pole_inputs = zeros(num_inputs, closed_pole_samples);
pole_outputs = zeros(num_outputs, closed_pole_samples);
pole_states(:, 1) = x0; %set the initial conditions

for k = 1:closed_pole_samples
    pole_inputs(:, k) = placed_controller * pole_states(:, k);              %u(k) = F * x(k) - linear control law
    pole_states(:, k + 1) = A * pole_states(:, k) + B * pole_inputs(:, k);  %x(k+1) = A*x(k) + B*u(k)
    pole_outputs(:, k) = C * pole_states(:, k) + D * pole_inputs(:, k);     %y(k) = C*x(k) + D*u(k)
end

if plot_all
    plot_two_mass('Pole Placement', 'Closed-Loop System with Pole Placement', pole_outputs, pole_inputs, update_file_path);   %function to plot positions and forces
end

%% Deadbeat Controller
deadbeat_poles = [-0.00001, 0.00001, 0.00000i, -0.00000i]; %plce will not allow for multiple poles at same location, so make all really close to 0
deadbeat_controller = -place(A, B, deadbeat_poles);
verify_poles = eig(A + B * deadbeat_controller);
if any(abs(verify_poles) > 1)
    printf('A Pole was placed outside of the unit circle, unstable controller!')
end

%Visualize the poles
if plot_all
    plot_pole_placement('Deadbeat Pole Placement', 'Deadbeat Pole Placement', verify_poles, deadbeat_poles, update_file_path);
end
deadbeat_samples = 20;
deadbeat_states = zeros(num_states, deadbeat_samples);
deadbeat_inputs = zeros(num_inputs, deadbeat_samples);
deadbeat_outputs = zeros(num_outputs, deadbeat_samples);
deadbeat_states(:, 1) = x0; %set the initial conditions

for k = 1:deadbeat_samples
    deadbeat_inputs(:, k) = deadbeat_controller * deadbeat_states(:, k);              %u(k) = F * x(k) - linear control law
    deadbeat_states(:, k + 1) = A * deadbeat_states(:, k) + B * deadbeat_inputs(:, k);  %x(k+1) = A*x(k) + B*u(k)
    deadbeat_outputs(:, k) = C * deadbeat_states(:, k) + D * deadbeat_inputs(:, k);     %y(k) = C*x(k) + D*u(k)
end

if plot_all
    plot_two_mass('Deadbeat Controller', 'Deadbeat Controller', deadbeat_outputs, deadbeat_inputs, update_file_path);
end

%% LQR Controller
Q = 100 * eye(num_states);  %cost of each states being away from 0
R = 1 * eye(num_inputs);    %cost od each input away from 0
gamma = 0.8;            %learning factor

%Consturct an LQR controller
%The LQR controller minimizes the cost function J = u' * R * u + x' * x
%Where u = F*x
F_lqr = discounted_LQR(A, B, gamma, Q, R);

lqr_samples = 200;
lqr_states = zeros(num_states, lqr_samples);
lqr_inputs = zeros(num_inputs, lqr_samples);
lqr_outputs = zeros(num_outputs, lqr_samples);
lqr_states(:, 1) = x0; %set the initial conditions

for k = 1:lqr_samples
    lqr_inputs(:, k) = F_lqr * lqr_states(:, k);              %u(k) = F * x(k) - linear control law
    lqr_states(:, k + 1) = A * lqr_states(:, k) + B * lqr_inputs(:, k);  %x(k+1) = A*x(k) + B*u(k)
    lqr_outputs(:, k) = C * lqr_states(:, k) + D * lqr_inputs(:, k);     %y(k) = C*x(k) + D*u(k)
end
%Find where it placed the poles
lqr_poles = eig(A + B * F_lqr);
if plot_all
    plot_pole_placement('Big Q LQR Pole Locations', sprintf('Q/R = %.d Pole Placement', Q(1, 1)/R(1, 1)), lqr_poles, -1, update_file_path);
end
%Plot full system out
if plot_all
    plot_two_mass('Big Q LQR Controller', sprintf('With LQR Controller (Q/R=%.d)', Q(1, 1)/R(1, 1)), lqr_outputs, lqr_inputs, update_file_path);
end

%% Demo Higher R LQR
R_big = 10 * eye(num_inputs);    %make relaive weightings more skewed to input
gamma = 0.8;


big_R_F_lqr = discounted_LQR(A, B, gamma, Q, R_big);

big_R_lqr_samples = 1200;
big_R_lqr_states = zeros(num_states, big_R_lqr_samples);
big_R_lqr_inputs = zeros(num_inputs, big_R_lqr_samples);
big_R_lqr_outputs = zeros(num_outputs, big_R_lqr_samples);
big_R_lqr_states(:, 1) = x0; %set the initial conditions

for k = 1:big_R_lqr_samples
    big_R_lqr_inputs(:, k) = big_R_F_lqr * big_R_lqr_states(:, k);              %u(k) = F * x(k) - linear control law
    big_R_lqr_states(:, k + 1) = A * big_R_lqr_states(:, k) + B * big_R_lqr_inputs(:, k);  %x(k+1) = A*x(k) + B*u(k)
    big_R_lqr_outputs(:, k) = C * big_R_lqr_states(:, k) + D * big_R_lqr_inputs(:, k);     %y(k) = C*x(k) + D*u(k)
end
%Find where it placed the poles
big_R_lqr_poles = eig(A + B * big_R_F_lqr);
if plot_all
    plot_pole_placement('Big R Pole Placement', sprintf('Q/R = %.d Pole Placement', Q(1, 1)/R_big(1, 1)), big_R_lqr_poles, -1, update_file_path);
end
%Plot full system out
if plot_all
    plot_two_mass('Big R LQR Controller', sprintf('With LQR Controller (Q/R = %.d)', Q(1, 1)/R_big(1, 1)), big_R_lqr_outputs, big_R_lqr_inputs, update_file_path);
end

%% ILC Introduction
%Setup ILC Parameters
p = 100; %number of steps to a process
num_ilc_states = num_outputs * p;  %effective states for ILC are the errors, generated from outputs
num_ilc_inputs = num_inputs * p;   %r inputs for every time step

x0_ilc = x0;

%Descriptive Matrix
[P_full, d_full] = P_from_ABCD(A, B, C, D, p, x0_ilc); %generate a descriptive matrix which captures the relation of inputs to outputs

%ILC State-space:
% e(j+1) = e(j) - P*del_u(j)
% y(j) = e(j)
A_ilc = eye(num_ilc_states);
B_ilc = -P_full;
C_ilc = eye(num_ilc_states);
D_ilc = 0;

if ~is_controllable(A_ilc, B_ilc) %always check it is even possible to control
    printf('The ILC System is not Controllable!')
end

F_ilc = 0.8 * pinv(P_full);     %Note: any controller that palces poles inside the unit circle works.

%Set the Goal Outputs (only have goals for p steps, regardless of lead)
%Use function draw_to_XY(p) if you would like to set your own
%drawing/special output. This is not written to be super robust, choose
%matchig p values
%[drawn_x, drawn_y] = draw_to_XY(p);
%save(sprintf('Saved Data\shape_p%.d.mat',p), 'drawn_x', 'drawn_y');    %save to a file so we only have to do this once

y_star_x = cos(2 * pi * (0:(p-1)) / p)';
y_star_y = sin(2 * pi * (0:(p-1)) / p)';
goal_matrix = [y_star_x, y_star_y]; %stack inputs next to eachother
y_star = reshape(goal_matrix', [], 1); %combine the seperate goals of each output into one vertical vector, alternating as necessary

%ILC Structure
num_trials = 10;    %number of trials to learn input
ILC_Trial(num_trials).output = [];
ILC_Trial(num_trials).input = [];
ILC_Trial(num_trials).del_u = [];   %relevant control parameter
ILC_Trial(num_trials).output_error = [];   %relevant control parameter

%First Trial
trial_num = 1;
ILC_Trial(trial_num).del_u = zeros(num_ilc_inputs, 1);  %first trial we set whatever del_u we want
ILC_Trial(trial_num).input = zeros(num_ilc_inputs, 1);  %Similarly, first input is arbitary

ILC_Trial(trial_num).output = [C * x0_ilc + D * ILC_Trial(trial_num).input(1:num_inputs); P_full * ILC_Trial(trial_num).input + d_full]; %simulate y(0) -> y(p + num_lead), since y = Pu + d ignores y(0), use IC
relevant_output = ILC_Trial(trial_num).output((num_outputs + 1):end);  %only compare to the ones we can / want to control
ILC_Trial(trial_num).output_error = y_star - relevant_output;

%Subsequent Iterations
for trial_num = 2:num_trials
    %Inputs
    ILC_Trial(trial_num).del_u = F_ilc * ILC_Trial(trial_num - 1).output_error;
    ILC_Trial(trial_num).input = ILC_Trial(trial_num).del_u + ILC_Trial(trial_num - 1).input;   %d_u(j) = u(j) - u(j-1)
    %Output
    ILC_Trial(trial_num).output = [C * x0_ilc + D * ILC_Trial(trial_num).input(1:num_inputs); P_full * ILC_Trial(trial_num).input + d_full];
    relevant_output = ILC_Trial(trial_num).output((1 + num_outputs):end);
    %Error
    ILC_Trial(trial_num).output_error = y_star - relevant_output;
end

if plot_all
    to_plot_ilc = [1, 2, 5, num_trials];
    plot_dual_ilc('Placed Controller ILC', 'ILC with a Placed Controller', ILC_Trial, y_star, -1, to_plot_ilc, update_file_path); %update_file_path
end

%% Arbitrary ILC
%Setup ILC Parameters
if false
    clear 'ILC_Trial' %make sure we do not carry over any ILC from previous
    p = 500; %number of steps to a process
    num_ilc_states = num_outputs * p;  %effective states for ILC are the errors, generated from outputs
    num_ilc_inputs = num_inputs * p;   %r inputs for every time step
    
    x0_ilc = x0;
    
    %Descriptive Matrix
    [P_full, d_full] = P_from_ABCD(A, B, C, D, p, x0_ilc); %generate a descriptive matrix which captures the relation of inputs to outputs
    
    A_ilc = eye(num_ilc_states);
    B_ilc = -P_full;
    C_ilc = eye(num_ilc_states);
    D_ilc = 0;
    
    if ~is_controllable(A_ilc, B_ilc) %always check it is even possible to control
        printf('The ILC System is not Controllable!')
    end
    
    F_ilc = 0.8 * pinv(P_full);     %Note: any controller that palces poles inside the unit circle works.
    
    loadedShape = load('Saved Data\dartmouth_p500.mat');    %read in the file
    y_star_x = loadedShape.drawn_x';
    y_star_y = loadedShape.drawn_y';
    scale = 10;
    goal_matrix = [y_star_x, y_star_y]; %stack inputs next to eachother
    y_star = scale * reshape(goal_matrix', [], 1); %combine the seperate goals of each output into one vertical vector, alternating as necessary
    
    %ILC Structure
    num_trials = 10;    %number of trials to learn input
    ILC_Trial(num_trials).output = [];
    ILC_Trial(num_trials).input = [];
    ILC_Trial(num_trials).del_u = [];   %relevant control parameter
    ILC_Trial(num_trials).output_error = [];   %relevant control parameter
    
    %First Trial
    trial_num = 1;
    ILC_Trial(trial_num).del_u = zeros(num_ilc_inputs, 1);  %first trial we set whatever del_u we want
    ILC_Trial(trial_num).input = zeros(num_ilc_inputs, 1);  %Similarly, first input is arbitary
    
    ILC_Trial(trial_num).output = [C * x0_ilc + D * ILC_Trial(trial_num).input(1:num_inputs); P_full * ILC_Trial(trial_num).input + d_full]; %simulate y(0) -> y(p + num_lead), since y = Pu + d ignores y(0), use IC
    relevant_output = ILC_Trial(trial_num).output((num_outputs + 1):end);  %only compare to the ones we can / want to control
    ILC_Trial(trial_num).output_error = y_star - relevant_output;
    
    %Subsequent Iterations
    for trial_num = 2:num_trials
        %Inputs
        ILC_Trial(trial_num).del_u = F_ilc * ILC_Trial(trial_num - 1).output_error;
        ILC_Trial(trial_num).input = ILC_Trial(trial_num).del_u + ILC_Trial(trial_num - 1).input;   %d_u(j) = u(j) - u(j-1)
        %Output
        ILC_Trial(trial_num).output = [C * x0_ilc + D * ILC_Trial(trial_num).input(1:num_inputs); P_full * ILC_Trial(trial_num).input + d_full];
        relevant_output = ILC_Trial(trial_num).output((1 + num_outputs):end);
        %Error
        ILC_Trial(trial_num).output_error = y_star - relevant_output;
    end
    
    if plot_all
        to_plot_ilc = [1, 2, 5, num_trials];
        
        resave_dartmouth = update_file_path
        %Super expliclity re-save the dartmouth plot. When doing so, go into
        %save_figure and adjust legend location for ONLY this one. Just
        %because where it draws things
    
        plot_dual_ilc('Dartmouth ILC', 'Arbitrary ILC with a Placed Controller', ILC_Trial, y_star, -1, to_plot_ilc, resave_dartmouth); %update_file_path
    end
else
    sprintf('Not doing the big ILC to save computation time')
end

%% Reinforcement Learning - Policy Iteration
converged_k = 20;   %number of steps to run out for with 'converged' controller
num_controllers = 5; %number of iterations of each controller before updating
Pj_dim = (num_states + num_inputs); %minimum dimensions of the P matrix such that we can solve (noise free)
num_collections = Pj_dim^2;  %number of data collections to make before solving the inv (noise free)

total_tries = num_controllers * num_collections + converged_k;

F_policy = zeros(num_inputs, num_states, 1);    %start with no controller

policy_states = zeros(num_states, total_tries);
policy_states(:, 1) = x0;
policy_inputs = zeros(num_inputs, total_tries);
policy_outputs = zeros(num_outputs, total_tries);
k = 1;


for iteration = 1:num_controllers
    %Reset Stacks each iteration
    Uk_stack = zeros(num_collections, 1);
    Xk_stack = zeros(num_collections, (Pj_dim)^2);

    %Compute stack infos
    for collection = 1:num_collections
        %Simulate nature
        %policy_inputs(:, k) = rand(num_inputs, 1);  %can replicate the 'random state' approach by just randomizing inputs
        policy_inputs(:, k) = F_policy(:, :, end) * policy_states(:, k) + rand_range(num_inputs, 1, -1, 1);    %compute input and exploration term

        %Next Step regardless of nature vs random
        policy_states(:, k + 1) = A * policy_states(:, k) + B * policy_inputs(:, k);
        policy_outputs(:, k) = C * policy_states(:, k) + D * policy_inputs(:, k);


        %Construct stacks
        xu_stack = [policy_states(:, k); policy_inputs(:, k)];
        xu_next_stack = [policy_states(:, k + 1); F_policy(:, :, end) * policy_states(:, k + 1)];
        Xk_stack(collection, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
        Uk_stack(collection, :) = policy_states(:, k)' * Q * policy_states(:, k) + policy_inputs(:, k)' * R * policy_inputs(:, k);  %where Q and R come into play

        %Move to next
        k = k + 1;
    end
    % [P, S, V] = svd(Xk_stack);
    % rank(S), iteration %we do not need the stack to be full rank, but it
    % must not drop ranks thorugh iterations
    %Calculate P and new controller
    PjS = pinv(Xk_stack) * Uk_stack;
    Pj = reshape(PjS, Pj_dim, Pj_dim);  %undo the stack
    Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
    Pjuu = Pj((num_states+1):end, (num_states+1):end);  %grab the bottom right quadrant, Puu
    PjxuT = Pj((num_states+1):end, 1:num_states);       %Grab bottom left, Pxu transpose
    new_F = -pinv(Pjuu) * PjxuT;    %definition of controller
    F_policy(:, :, end + 1) = new_F;    %add to end of list
end

%Plot / simulate past convergance (without exploration/noise) so actual
%control possible
for ndx = 1:converged_k
    %Simulate nature
    policy_inputs(:, k) = F_policy(:, :, end) * policy_states(:, k);
    policy_states(:, k + 1) = A * policy_states(:, k) + B * policy_inputs(:, k);
    policy_outputs(:, k) = C * policy_states(:, k) + D * policy_inputs(:, k);
    k = k + 1;
end

if plot_all
    plot_two_mass('Policy Iteration IO', 'Under Policy Iteration', policy_outputs, policy_inputs, update_file_path);
    plot_controller_history('Policy Iteration Controller', 'Under Policy Iteration', F_policy, -1, -1, update_file_path);   %-1s denote to plot every state and input weight
end


%% Reinforcement Learning - Input Decoupling
converged_k = 20;
num_controllers = 5;
Pj_dim = (num_states + 1); %effective input count is 1 now
num_collections = Pj_dim^2;

total_tries = num_controllers * num_collections + converged_k;

F_decoupled = zeros(num_inputs, num_states, 1);    %start with no controller

decoupled_states = zeros(num_states, total_tries);
decoupled_states(:, 1) = x0;
decoupled_inputs = zeros(num_inputs, total_tries);
decoupled_outputs = zeros(num_outputs, total_tries);
k = 1;

for iteration = 1:num_controllers
    for input_num = 1:num_inputs
        %Reset Stacks each iteration
        Uk_stack = zeros(num_collections, 1);
        Xk_stack = zeros(num_collections, (Pj_dim)^2);
        F_i = F_decoupled(input_num, :, end);   %input on the current input of inspection
        %Compute stack infos
        for collection = 1:num_collections
            %Simulate nature
            %decoupled_inputs(:, k) = rand(num_inputs, 1);  %can replicate the 'random state' approach by just randomizing inputs
            decoupled_inputs(:, k) = F_decoupled(:, :, end) * decoupled_states(:, k);    %compute input
            decoupled_inputs(input_num, k) = rand_range(1, 1, -1, 1); %+decoupled_inputs(input_num, k); %

            %Next Step regardless of nature vs random
            decoupled_states(:, k + 1) = A * decoupled_states(:, k) + B * decoupled_inputs(:, k);
            decoupled_outputs(:, k) = C * decoupled_states(:, k) + D * decoupled_inputs(:, k);


            %Construct stacks
            xu_stack = [decoupled_states(:, k); decoupled_inputs(input_num, k)];
            xu_next_stack = [decoupled_states(:, k + 1); F_i * decoupled_states(:, k + 1)];
            Xk_stack(collection, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
            Uk_stack(collection, :) = decoupled_states(:, k)' * Q * decoupled_states(:, k) + decoupled_inputs(:, k)' * R * decoupled_inputs(:, k);  %where Q and R come into play

            %Move to next
            k = k + 1;
        end

        %Calculate P and new controller
        PjS = pinv(Xk_stack) * Uk_stack;
        Pj = reshape(PjS, Pj_dim, Pj_dim);  %undo the stack
        Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
        Pjuu = Pj((num_states+1):end, (num_states+1):end);  %grab the bottom right quadrant, Puu
        PjxuT = Pj((num_states+1):end, 1:num_states);       %Grab bottom left, Pxu transpose
        new_F_i = -pinv(Pjuu) * PjxuT;    %definition of controller
        F_decoupled(:, :, end + 1) = F_decoupled(:, :, end);    %copy old controller
        F_decoupled(input_num, :, end) = new_F_i;    %update current input
    end
end

%Plot / simulate past convergance (without exploration/noise) so actual
%control possible
for ndx = 1:converged_k
    %Simulate nature
    decoupled_inputs(:, k) = F_decoupled(:, :, end) * decoupled_states(:, k);
    decoupled_states(:, k + 1) = A * decoupled_states(:, k) + B * decoupled_inputs(:, k);
    decoupled_outputs(:, k) = C * decoupled_states(:, k) + D * decoupled_inputs(:, k);
    k = k + 1;
end

if plot_all
    plot_two_mass('Input Decoupling IO', 'Under Input Decoupling', decoupled_outputs, decoupled_inputs, update_file_path);
    plot_controller_history('Input Decoupling Controller', 'Under Input Decoupling', F_decoupled, -1, -1, update_file_path, true);   %note only updates every other controller #
end


%% Compute Costs of Different Controllers
cost_max_k = 500;
%LQR
cost_lqr_states = zeros(num_states, cost_max_k);
cost_lqr_inputs = zeros(num_inputs, cost_max_k);
%Policy
cost_policy_states = zeros(num_states, cost_max_k);
cost_policy_inputs = zeros(num_inputs, cost_max_k);
%Decoupled
cost_decoupled_states = zeros(num_states, cost_max_k);
cost_decoupled_inputs = zeros(num_inputs, cost_max_k);

%ICs
cost_decoupled_states(:, 1) = x0;
cost_policy_states(:, 1) = x0;
cost_lqr_states(:, 1) = x0;

%Cost values
lqr_cost = 0;
policy_cost = 0;
decoupled_cost = 0;

%Simulate
for k = 1:cost_max_k
    %LQR
    cost_lqr_inputs(:, k) = F_lqr * cost_lqr_states(:, k);
    cost_lqr_states(:, k + 1) = A * cost_lqr_states(:, k) + B * cost_lqr_inputs(:, k);
    lqr_cost = lqr_cost + (cost_lqr_inputs(:, k)' * R * cost_lqr_inputs(:, k) + cost_lqr_states(:, k)' * Q * cost_lqr_states(:, k));
    %Policy
    cost_policy_inputs(:, k) = F_policy(:, :, end) * cost_policy_states(:, k);
    cost_policy_states(:, k + 1) = A * cost_policy_states(:, k) + B * cost_policy_inputs(:, k);
    policy_cost = policy_cost + (cost_policy_inputs(:, k)' * R * cost_policy_inputs(:, k) + cost_policy_states(:, k)' * Q * cost_policy_states(:, k));
    %Decoupled
    cost_decoupled_inputs(:, k) = F_decoupled(:, :, end) * cost_decoupled_states(:, k);
    cost_decoupled_states(:, k + 1) = A * cost_decoupled_states(:, k) + B * cost_decoupled_inputs(:, k);
    decoupled_cost = decoupled_cost + (cost_decoupled_inputs(:, k)' * R * cost_decoupled_inputs(:, k) + cost_decoupled_states(:, k)' * Q * cost_decoupled_states(:, k));
end

%% Compare LQR, Policy, and Input Decoupled Controller
lqr_cost, policy_cost, decoupled_cost %all controllers should have same cost
F_lqr
F_policy(:, :, end)
F_decoupled(:, :, end)
norm(F_lqr - F_policy(:, :, end))/numel(F_lqr)
norm(F_lqr - F_decoupled(:, :, end))/numel(F_lqr)
