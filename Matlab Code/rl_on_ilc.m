%Appllication of Reinforcement Learning to the ILC Problem
%Noah Dunleavy
%Honors Thesis under direction of Professor Minh Q. Phan
%Thayer School of Engineering
clc; clear;
addpath('Saved Data', 'Functions');
setDefaultFigProp();
plot_policy = false;   %majority of run time goes to figure gen, so this speeds up when false
plot_decoupled = false;
update_file_path = -1;%'C:\Users\noahd\OneDrive\Desktop\Thesis\Thesis Images\RL on ILC';  %set this to the save location if want to update figures, or -1 if not
if update_file_path ~= -1
    keyboard    %ensure we have to be very concious about ever updating all the images
end

%% System Creation
thesis_system = load('Saved Data\thesis_system.mat');
A = thesis_system.A;
B = thesis_system.B;
C = thesis_system.C;
D = thesis_system.D;
num_inputs = thesis_system.num_inputs;
num_outputs = thesis_system.num_outputs;
num_states = thesis_system.num_states;
x0 = thesis_system.x0;

%% Goal Definition
p = 10;

% loadedShape = load('Saved Data\heart_p20.mat');    %read in the file
% y_star_x = loadedShape.drawn_x';
% y_star_y = loadedShape.drawn_y';

y_star_x = cos(2 * pi * (0:(p-1)) / p)';
y_star_y = sin(2 * pi * (0:(p-1)) / p)';
goal_matrix = [y_star_x, y_star_y]; %stack inputs next to eachother
y_star = reshape(goal_matrix', [], 1); %combine the seperate goals of each output into one vertical vector, alternating as necessary

%% ILC System
[P, d] = P_from_ABCD(A, B, C, D, p, x0);    %y(1:p) = P * u(0:(p-1)) + d
num_ilc_states = p * num_outputs;   %error bar - one for each output at each time step
num_ilc_inputs = p * num_inputs;    %u bar - each time step gets an input

%ILC State-space so that e(j+1) = A * e(j) + B * del_u_j1 and y(j) = e(j)
A_ilc = eye(num_ilc_states);
B_ilc = -P;
C_ilc = eye(num_ilc_states);
D_ilc = 0;
if (~is_controllable(A_ilc, B_ilc))
    fprintf('The ILC System is not Controllable!')
    return
end

%% Policy Iteration RL Parameters, R = 1
%Learning Weights
Q = 100 * eye(num_ilc_states);  %cost of each error
R = 1 * eye(num_ilc_inputs);    %penalize inputs, or more accurately change in input
    %Setting this too small causes controller to fail, but logically 0
    %should be safe
gamma = 0.8;
policy_exploration_mag = 1;
F_ilc_lqr = discounted_LQR(A_ilc, B_ilc, gamma, Q, R);  %perfect nowledge controller



%Iteration Counts
Pj_dim = (num_ilc_inputs + num_ilc_states); %the square diension of the Pj matrix
num_controllers = 5;    %number of controllers to try to create
num_collections_per_controller = Pj_dim^2;    %number of datasets needed per controller to update
num_converged = 1000;    %how many trials to then go through once we are done 'learning' (so we arent noisey)

%To be able to set R smaller, increase exploration magnitude ad the number
%of collections per controller

total_trial_count = num_controllers * num_collections_per_controller + num_converged;

%Preallocate the space and define structure
ILC_Trial(total_trial_count).output = [];   
ILC_Trial(total_trial_count).input = [];
ILC_Trial(total_trial_count).output_error = [];
ILC_Trial(total_trial_count).del_u = [];

%% Policy Iteration Learning Process, R = 1
rng('default'); rng(10);%ensure same randomization each time
F_ilc = zeros(num_ilc_inputs, num_ilc_states, num_controllers + 1);   %start with no controller

%Prepopulate the first trial
trial_num = 1;
ILC_Trial(trial_num).input = zeros(num_ilc_inputs, 1);  %start with open-loop / no input
ILC_Trial(trial_num).output = [C*x0; d];    %open loop response is IC and then d term
ILC_Trial(trial_num).output_error = y_star - d;    %relevant error

trial_num = 2; %start at second trial now
for iteration = 1:num_controllers    
    Uk_stack = zeros(num_collections_per_controller, 1);
    Xk_stack = zeros(num_collections_per_controller, (Pj_dim)^2);

    %Simulate the necessary trials
    for trial = 1:num_collections_per_controller   %number of trials to collect before updatin controller
        %ILC / Real Controller Process
        %Inputs
        exploration_term = rand_range(num_ilc_inputs, 1, -policy_exploration_mag, policy_exploration_mag);    %jiggle to learn
        ILC_Trial(trial_num).del_u = F_ilc(:, :, iteration) * ILC_Trial(trial_num - 1).output_error + exploration_term;
        ILC_Trial(trial_num).input = ILC_Trial(trial_num - 1).input + ILC_Trial(trial_num).del_u;
        %Simulate Reality
        relevant_output = P * ILC_Trial(trial_num).input + d;   %y(1) -> y(p)
        ILC_Trial(trial_num).output = [C*x0 + D*ILC_Trial(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
        %Calculate Error
        ILC_Trial(trial_num).output_error = y_star - relevant_output;
        %error_next_law = ILC_Trial(trial_num - 1).output_error - P * ILC_Trial(trial_num).del_u;   %verify that this matches the produced error, that is:
            %e(j+1) = e(j) - P*del_u(j+1) when del_u(j+1) = L * e(j) = u(j) - u(j-1)

        %RL Translation
        state = ILC_Trial(trial_num - 1).output_error; %analogous state, x(k) -> e_(j-1)
        input = ILC_Trial(trial_num).del_u;%analogous input, u(k) -> del_u_(j) = L * e_(j-1)
        next_state = ILC_Trial(trial_num).output_error;    %x(k+1) = e(j)
        next_input = F_ilc(:, :, iteration) * next_state; %no exploration term here

        xu_stack = [state; input];
        xu_next_stack = [next_state; next_input];

        Xk_stack(trial, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
        Uk_stack(trial, :) = input' * R * input + state' * Q * state;

        trial_num = trial_num + 1;
    end

    %Calculate P and new controller
    % [U, S, V] = svd(Xk_stack);
    % rank(S)
    PjS = pinv(Xk_stack) * Uk_stack;
    Pj = reshape(PjS, Pj_dim, Pj_dim);
    Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
    Pjuu = Pj((num_ilc_states+1):end, (num_ilc_states+1):end);
    PjxuT = Pj((num_ilc_states+1):end, 1:num_ilc_states);
    new_F = -pinv(Pjuu) * PjxuT;
    F_ilc(:, :, iteration + 1) = new_F;
end

for ndx = 1:num_converged
    %Inputs
    ILC_Trial(trial_num).del_u = F_ilc(:, :, end) * ILC_Trial(trial_num - 1).output_error;
    ILC_Trial(trial_num).input = ILC_Trial(trial_num - 1).input + ILC_Trial(trial_num).del_u;
    %Simulate Reality
    relevant_output = P * ILC_Trial(trial_num).input + d;   %y(1) -> y(p)
    ILC_Trial(trial_num).output = [C*x0 + D*ILC_Trial(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
    %Calculate Error
    ILC_Trial(trial_num).output_error = y_star - relevant_output;
    trial_num = trial_num + 1;
end
controller_error = norm(F_ilc(:, :, end) - F_ilc_lqr)/numel(F_ilc_lqr)

%% Visualize Learning
if plot_policy
    to_plot = 4;
    plot_dual_ilc('Policy Iteration on ILC', 'Policy Iteration for ILC', ILC_Trial, y_star, -1, to_plot, update_file_path);
    plot_controller_history('Policy Iteration ILC Controller', 'Policy Iteration ILC Controller Weights', F_ilc, 4, 5, update_file_path); %we cannot plot *all* the IOs because there are so many, pick a few
end

%% Input Decoupled Learning
%Iteration Counts
decoupled_exploration_mag = policy_exploration_mag;
Pj_dim = (1 + num_ilc_states); 
num_controllers = 5;    
num_collections_per_controller = Pj_dim^2;   
num_converged_decoupled = num_collections_per_controller;    

total_trial_count_decoupled = num_ilc_inputs * num_controllers * num_collections_per_controller + num_converged_decoupled; %now need a * num inputs

%Preallocate the space and define structure
ILC_Trial_decoupled(total_trial_count_decoupled).output = [];   
ILC_Trial_decoupled(total_trial_count_decoupled).input = [];
ILC_Trial_decoupled(total_trial_count_decoupled).output_error = [];
ILC_Trial_decoupled(total_trial_count_decoupled).del_u = [];

%% Input Decoupled Learning Process
rng('default'); rng(10);%ensure same randomization each time
F_ilc_decoupled = zeros(num_ilc_inputs, num_ilc_states, 1);   

%Prepopulate the first trial
trial_num = 1;
ILC_Trial_decoupled(trial_num).input = zeros(num_ilc_inputs, 1);  %start with open-loop / no input
ILC_Trial_decoupled(trial_num).output = [C*x0; d];    %open loop response is IC and then d term
ILC_Trial_decoupled(trial_num).output_error = y_star - d;    %relevant error

trial_num = 2; %start at second trial now
for iteration = 1:(num_controllers)    
    for input_num = 1:num_ilc_inputs
        Uk_stack = zeros(num_collections_per_controller, 1);
        Xk_stack = zeros(num_collections_per_controller, (Pj_dim)^2);
        F_i = F_ilc_decoupled(input_num, :, end);
        %Simulate the necessary trials
        for trial = 1:num_collections_per_controller   %number of trials to collect before updatin controller
            %ILC / Real Controller Process
            %Inputs
            exploration_term = rand_range(1, 1, -decoupled_exploration_mag, decoupled_exploration_mag);    %jiggle to learn
            ILC_Trial_decoupled(trial_num).del_u = F_ilc_decoupled(:, :, iteration) * ILC_Trial_decoupled(trial_num - 1).output_error;

            ILC_Trial_decoupled(trial_num).del_u(input_num) = ILC_Trial_decoupled(trial_num).del_u(input_num) + exploration_term; %differnece between PI - only learn on one input

            ILC_Trial_decoupled(trial_num).input = ILC_Trial_decoupled(trial_num - 1).input + ILC_Trial_decoupled(trial_num).del_u;
            %Simulate Reality
            relevant_output = P * ILC_Trial_decoupled(trial_num).input + d;   %y(1) -> y(p)
            ILC_Trial_decoupled(trial_num).output = [C*x0 + D*ILC_Trial_decoupled(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
            %Calculate Error
            ILC_Trial_decoupled(trial_num).output_error = y_star - relevant_output;


            %RL Translation
            state = ILC_Trial_decoupled(trial_num - 1).output_error; %analogous state, x(k) -> e_(j-1)
            full_input = ILC_Trial_decoupled(trial_num).del_u;  %Uk vs Xk use different inputs
            input = full_input(input_num);
            next_state = ILC_Trial_decoupled(trial_num).output_error;    %x(k+1) = e(j)
            next_input = F_i * next_state; %no exploration term here

            xu_stack = [state; input];
            xu_next_stack = [next_state; next_input];

            Xk_stack(trial, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
            Uk_stack(trial, :) = full_input' * R * full_input + state' * Q * state;

            trial_num = trial_num + 1;
        end

        %Calculate P and new controller
        PjS = pinv(Xk_stack) * Uk_stack;
        Pj = reshape(PjS, Pj_dim, Pj_dim);
        Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
        Pjuu = Pj((num_ilc_states+1):end, (num_ilc_states+1):end);
        PjxuT = Pj((num_ilc_states+1):end, 1:num_ilc_states);
        new_F_i = -pinv(Pjuu) * PjxuT;
        F_ilc_decoupled(:, :, end + 1) = F_ilc_decoupled(:, :, end);
        F_ilc_decoupled(input_num, :, end) = new_F_i;
    end
end

for ndx = 1:num_converged_decoupled
    %Inputs
    ILC_Trial_decoupled(trial_num).del_u = F_ilc_decoupled(:, :, end) * ILC_Trial_decoupled(trial_num - 1).output_error;
    ILC_Trial_decoupled(trial_num).input = ILC_Trial_decoupled(trial_num - 1).input + ILC_Trial_decoupled(trial_num).del_u;
    %Simulate Reality
    relevant_output = P * ILC_Trial_decoupled(trial_num).input + d;   %y(1) -> y(p)
    ILC_Trial_decoupled(trial_num).output = [C*x0 + D*ILC_Trial_decoupled(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
    %Calculate Error
    ILC_Trial_decoupled(trial_num).output_error = y_star - relevant_output;
    trial_num = trial_num + 1;
end
controller_error = norm(F_ilc_decoupled(:, :, end) - F_ilc_lqr)/numel(F_ilc_lqr)

%% Visualize Learning
if plot_decoupled
    to_plot = 4;
    plot_dual_ilc('Input Decoupling on ILC', 'Input Decoupling for ILC', ILC_Trial_decoupled, y_star, -1, to_plot, update_file_path);
    plot_controller_history('Input Decoupling ILC Controller', 'Input Decoupling ILC Controller Weights', F_ilc_decoupled, 4, 5, update_file_path, true);
end

%% Policy Iteration RL Parameters, R = 1e-6
%Learning Weights
Q = 100 * eye(num_ilc_states);  %cost of each error
small_R = 1e-6 * eye(num_ilc_inputs);    %penalize inputs, or more accurately change in input
    %Setting this too small causes controller to fail, but logically 0
    %should be safe
gamma = 0.8;
small_R_policy_exploration_mag = 1000;
F_ilc_lqr_small_R = discounted_LQR(A_ilc, B_ilc, gamma, Q, small_R);  %perfect nowledge controller

%Iteration Counts
Pj_dim = (num_ilc_inputs + num_ilc_states); %the square diension of the Pj matrix
num_controllers = 5;    %number of controllers to try to create
num_collections_per_controller = Pj_dim^2;    %number of datasets needed per controller to update
num_converged = 100;    %how many trials to then go through once we are done 'learning' (so we arent noisey)

%To be able to set R smaller, increase exploration magnitude ad the number
%of collections per controller

total_trial_count = num_controllers * num_collections_per_controller + num_converged;

%Preallocate the space and define structure
ILC_Trial_small_R(total_trial_count).output = [];   
ILC_Trial_small_R(total_trial_count).input = [];
ILC_Trial_small_R(total_trial_count).output_error = [];
ILC_Trial_small_R(total_trial_count).del_u = [];

%% Policy Iteration Learning Process, R = 1e-6
rng('default'); rng(10);%ensure same randomization each time
F_ilc_policy_small_R = zeros(num_ilc_inputs, num_ilc_states, num_controllers + 1);   %start with no controller

%Prepopulate the first trial
trial_num = 1;
ILC_Trial_small_R(trial_num).input = zeros(num_ilc_inputs, 1);  %start with open-loop / no input
ILC_Trial_small_R(trial_num).output = [C*x0; d];    %open loop response is IC and then d term
ILC_Trial_small_R(trial_num).output_error = y_star - d;    %relevant error

trial_num = 2; %start at second trial now
for iteration = 1:num_controllers    
    Uk_stack = zeros(num_collections_per_controller, 1);
    Xk_stack = zeros(num_collections_per_controller, (Pj_dim)^2);

    %Simulate the necessary trials
    for trial = 1:num_collections_per_controller   %number of trials to collect before updatin controller
        %ILC / Real Controller Process
        %Inputs
        exploration_term = rand_range(num_ilc_inputs, 1, -small_R_policy_exploration_mag, small_R_policy_exploration_mag);    %jiggle to learn
        ILC_Trial_small_R(trial_num).del_u = F_ilc_policy_small_R(:, :, iteration) * ILC_Trial_small_R(trial_num - 1).output_error + exploration_term;
        ILC_Trial_small_R(trial_num).input = ILC_Trial_small_R(trial_num - 1).input + ILC_Trial_small_R(trial_num).del_u;
        %Simulate Reality
        relevant_output = P * ILC_Trial_small_R(trial_num).input + d;   %y(1) -> y(p)
        ILC_Trial_small_R(trial_num).output = [C*x0 + D*ILC_Trial_small_R(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
        %Calculate Error
        ILC_Trial_small_R(trial_num).output_error = y_star - relevant_output;
        %error_next_law = ILC_Trial(trial_num - 1).output_error - P * ILC_Trial(trial_num).del_u;   %verify that this matches the produced error, that is:
            %e(j+1) = e(j) - P*del_u(j+1) when del_u(j+1) = L * e(j) = u(j) - u(j-1)

        %RL Translation
        state = ILC_Trial_small_R(trial_num - 1).output_error; %analogous state, x(k) -> e_(j-1)
        input = ILC_Trial_small_R(trial_num).del_u;%analogous input, u(k) -> del_u_(j) = L * e_(j-1)
        next_state = ILC_Trial_small_R(trial_num).output_error;    %x(k+1) = e(j)
        next_input = F_ilc_policy_small_R(:, :, iteration) * next_state; %no exploration term here

        xu_stack = [state; input];
        xu_next_stack = [next_state; next_input];

        Xk_stack(trial, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
        Uk_stack(trial, :) = input' * small_R * input + state' * Q * state;

        trial_num = trial_num + 1;
    end

    %Calculate P and new controller
    % [U, S, V] = svd(Xk_stack);
    % rank(S)
    PjS = pinv(Xk_stack) * Uk_stack;
    Pj = reshape(PjS, Pj_dim, Pj_dim);
    Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
    Pjuu = Pj((num_ilc_states+1):end, (num_ilc_states+1):end);
    PjxuT = Pj((num_ilc_states+1):end, 1:num_ilc_states);
    new_F = -pinv(Pjuu) * PjxuT;
    F_ilc_policy_small_R(:, :, iteration + 1) = new_F;
end

for ndx = 1:num_converged
    %Inputs
    ILC_Trial_small_R(trial_num).del_u = F_ilc_policy_small_R(:, :, end) * ILC_Trial_small_R(trial_num - 1).output_error;
    ILC_Trial_small_R(trial_num).input = ILC_Trial_small_R(trial_num - 1).input + ILC_Trial_small_R(trial_num).del_u;
    %Simulate Reality
    relevant_output = P * ILC_Trial_small_R(trial_num).input + d;   %y(1) -> y(p)
    ILC_Trial_small_R(trial_num).output = [C*x0 + D*ILC_Trial_small_R(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
    %Calculate Error
    ILC_Trial_small_R(trial_num).output_error = y_star - relevant_output;
    trial_num = trial_num + 1;
end
controller_error = norm(F_ilc_policy_small_R(:, :, end) - F_ilc_lqr_small_R)/numel(F_ilc_lqr_small_R)

%% Visualize Learning
if plot_policy
    to_plot = 4;
    plot_dual_ilc('Small R Policy Iteration on ILC', 'Policy Iteration for ILC - Reduced R', ILC_Trial_small_R, y_star, -1, to_plot, update_file_path);
    plot_controller_history('Small R Policy Iteration ILC Controller', 'Policy Iteration ILC Controller Weights - Reduced R', F_ilc_policy_small_R, 4, 5, update_file_path); %we cannot plot *all* the IOs because there are so many, pick a few
end
