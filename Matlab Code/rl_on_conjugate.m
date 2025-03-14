%Appllication of RL to the ILC Problem with Conjugate
%Basis Functions
%Noah Dunleavy
%Honors Thesis under direction of Professor Minh Q. Phan
%Thayer School of Engineering
clc; clear;
addpath('Saved Data', 'Functions');
setDefaultFigProp();
plot_all = false;
update_file_path = -1;%'C:\Users\noahd\OneDrive\Desktop\Thesis\Thesis Images\RL on Conjugate';  %set this to the save location if want to update figures, or -1 if not
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

%% ILC Parameters
p = 10;

num_ilc_states = p * num_outputs;
num_ilc_inputs = p * num_inputs;

[P, d] = P_from_ABCD(A, B, C, D, p, x0);

%% Goal Definition
y_star_x = cos(2 * pi * (0:(p-1)) / p)';
y_star_y = sin(2 * pi * (0:(p-1)) / p)';
goal_matrix = [y_star_x, y_star_y]; %stack inputs next to eachother
y_star = reshape(goal_matrix', [], 1); %combine the seperate goals of each output into one vertical vector, alternating as necessary

%% Conjugate Basis Creation (Batch)
%Learning Weights for the Batch Learning
Q_batch = 100*eye(num_ilc_states);
R_batch = 0 * eye(num_ilc_inputs);
max_conj = num_ilc_inputs;
[conjugate_basis_functions, conjugate_betas] = generate_conjugate(num_ilc_inputs, max_conj, P, Q_batch, R_batch, d, y_star);

%% Conjugate Basis Functions Definition
num_basis_output = max_conj; %full definition / resolution
num_basis_input = num_basis_output;

output_basis_functions = conjugate_basis_functions; %full resolution output %conjugate_basis_functions;
input_basis_functions = conjugate_basis_functions;

output_basis_functions_pinv = pinv(output_basis_functions);   %since done a lot, just compute once and use as constant

%% Render Goal
if (plot_all)
    fig = figure('Name', 'Goal Output for Arbitrary Basis ILC');
    plot(y_star(1:2:end), y_star(2:2:(end+1)))
    axis equal
    %Should have already saved in RL on ILC
    %save_figure(update_file_path, fig);
end

%% F_lqr with Varying Basis Functions
%RL Parameters
Q = 100 * eye(num_ilc_states);
R = 10 * eye(num_ilc_inputs);
gamma = 0.8;


output_basis = output_basis_functions;
A_ilc_basis = eye(num_ilc_states); %fixed - capture the whole output
for ndx = 1:3
    current_basis = ndx;
    if ndx == 3
        current_basis = [1, 2];
    end
    input_basis = conjugate_basis_functions(:, current_basis);
    H = -pinv(output_basis) * P * input_basis;
    F_lqr = discounted_LQR(A_ilc_basis, -H, gamma, Q, R(current_basis, current_basis));
end

%% Full Resolution Demo
FIFO_Trials = basis_ilc_sim(P, d, C*x0, conjugate_basis_functions, conjugate_basis_functions, y_star, 20);
if plot_all
    plot_dual_ilc('FIFO Conjugate Basis', 'Full Resolution Conjugate Basis', FIFO_Trials, y_star, -1, 4, update_file_path);
end

%% Fixed Resolution T_y = I
%Fix the output basis
num_output_basis = num_ilc_inputs; %assume worst case
output_basis_functions = eye(num_output_basis);
%Conjugate episodes holder
Episode = []; 
%Setup learning parameters
Q = 100 * eye(num_output_basis);
big_R = 1 * eye(num_ilc_inputs);
small_R = 1e-6 * eye(num_ilc_inputs);
gamma = 0.8;
exploration_mag = 10;
num_controllers = 5;
num_converged = 10;
fixed_I_Trials_big_R = []; %hold all the trials
fixed_I_Trials_small_R = []; 
F_big_R = zeros(num_ilc_inputs, num_output_basis);
F_small_R = zeros(num_ilc_inputs, num_output_basis);
for phi_num = 1:num_ilc_inputs %try every combo possible
    %Update basis space
    Episode = generate_iterative_conjugate(P, d, y_star, Episode, Q, 0);
    phi = Episode(phi_num).Phi(:, phi_num); %get the most recent phi
    %Learn the controller with big R
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q, big_R(phi_num, phi_num), num_controllers, exploration_mag, phi, output_basis_functions, num_converged);
    F_big_R(phi_num, :) = F_phi(:, :, end); %save the controller in the stack
    fixed_I_Trials_big_R = [fixed_I_Trials_big_R, temp_Trial]; %save teh trial info
    
    %Controller with small R
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q, small_R(phi_num, phi_num), num_controllers, exploration_mag, phi, output_basis_functions, num_converged);
    F_small_R(phi_num, :) = F_phi(:, :, end); %save the controller in the stack
    fixed_I_Trials_small_R = [fixed_I_Trials_small_R, temp_Trial]; %save teh trial info
end
F_lqr_big_R = discounted_LQR(eye(num_output_basis), -pinv(output_basis_functions) * P * Episode(end).Phi, gamma, Q, big_R);
F_lqr_small_R = discounted_LQR(eye(num_output_basis), -pinv(output_basis_functions) * P * Episode(end).Phi, gamma, Q, small_R);


%Show ow they compare
sprintf('F_lqr vs Policy when Fixed Output Basis of I with Big R')
lqr = F_lqr_big_R([1, 2, 19, 20], [1, 2, 19, 20])
policy = F_big_R([1, 2, 19, 20], [1, 2, 19, 20])

sprintf('F_lqr vs Policy when Fixed Output Basis of I with Small R')
lqr = F_lqr_small_R([1, 2, 19, 20], [1, 2, 19, 20])
policy = F_small_R([1, 2, 19, 20], [1, 2, 19, 20])

%Plot the path through learning and straight application
if plot_all
    plot_dual_ilc('Conjugate Phi with Fixed I Output Basis', 'Learning with Singular Input Basis', fixed_I_Trials_big_R, y_star, -1, 4, update_file_path);

    one_trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi(:, 1), output_basis_functions, y_star, 10, F_big_R(1, :)); 
    plot_dual_ilc('Controller Application of One Conjugate Input Basis with Fixed I Output Basis', 'One-Controller Fixed I Output Basis', one_trial, y_star, -1, 4, update_file_path);


    half_trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi(:, 1:10), output_basis_functions, y_star, 10, F_big_R(1:10, :)); 
    plot_dual_ilc('Half Controller Application of Conjugate Input Basis with Fixed I Output Basis', 'Half-Controller Fixed I Output Basis', half_trial, y_star, -1, 4, update_file_path);

    seven_five_trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi(:, 1:15), output_basis_functions, y_star, 10, F_big_R(1:15, :)); 
    plot_dual_ilc('Three-Quarter Controller Application of Conjugate Input Basis with Fixed I Output Basis', 'Three-Quarter Fixed I Output Basis', seven_five_trial, y_star, -1, 4, update_file_path);

    final_Trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi, output_basis_functions, y_star, 10, F_big_R); 
    plot_dual_ilc('Controller Application of Conjugate Input Basis with Fixed I Output Basis', 'Fixed I Output Basis', final_Trial, y_star, -1, 4, update_file_path);
end

%% Fixed Resolution T_y = y*
%Fix the output basis
num_output_basis = 1; 
output_basis_functions = 100 * y_star;
%Conjugate episodes holder
clear 'Episode';
Episode = []; 
%Setup learning parameters
Q = 1000 * eye(num_output_basis);
R_y_star = 10 * eye(num_ilc_inputs);
gamma = 0.8;
exploration_mag = 1;
num_controllers = 2;
num_converged = 10;
y_star_fixed_I_Trials = []; %hold all the trials
F = zeros(num_ilc_inputs, num_output_basis);
for phi_num = 1:num_ilc_inputs %try every combo possible
    %Update basis space
    Episode = generate_iterative_conjugate(P, d, y_star, Episode, Q, 0);
    phi = Episode(phi_num).Phi(:, phi_num); %get the most recent phi
    %Learn the controller
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q, R_y_star(phi_num, phi_num), num_controllers, exploration_mag, phi, output_basis_functions, num_converged);
    F(phi_num, :) = F_phi(:, :, end); %save the controller in the stack
    y_star_fixed_I_Trials = [y_star_fixed_I_Trials, temp_Trial]; %save teh trial info
end
F_lqr_y_star = discounted_LQR(eye(num_output_basis), -pinv(output_basis_functions) * P * Episode(end).Phi, gamma, Q, R_y_star);


%Show ow they compare
sprintf('F_lqr vs Policy when Fixed Output Basis of y*')
lqr = F_lqr_y_star([1, 2, 19, 20])
policy = F([1, 2, 19, 20])

num_sim = 100;
y_star_Trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi, output_basis_functions, y_star, num_sim, F); 
F_lqr_y_star_mod = discounted_LQR(eye(num_output_basis), -pinv(output_basis_functions) * P * Episode(end).Phi, gamma, Q, R_y_star*(1e-8));
y_star_lqr_Trial = basis_ilc_sim(P, d, C*x0, Episode(end).Phi, output_basis_functions, y_star, num_sim, F_lqr_y_star_mod); 
if plot_all
    plot_dual_ilc('Controller Application of Conjugate Input Basis with y star Output Basis', 'Fixed y* Output Basis RL', y_star_Trial, y_star, -1, 4, update_file_path);
    plot_dual_ilc('LQR Controller Application of Conjugate Input Basis with y star Output Basis', 'Fixed y* Output Basis LQR', y_star_lqr_Trial, y_star, -1, 4, update_file_path);
end

%% Fixed Resolution T_y = Phi^b

fixed_phi_R = eye(num_ilc_inputs);
exploration_mag = 1; %little extra noise to be safe
T_y_scale = 10;

num_converged = 0;
num_controllers = 5;
%Generate our conjugate basis
conjugate_basis_functions = Episode(end).Phi; %use our previous iteratively calculated instead of doing a batch - computationally faster and more consistent
%Fix the output basis
num_output_basis = num_ilc_inputs; %assume worst case
output_basis_functions = conjugate_basis_functions(:, 1:num_output_basis) * T_y_scale;
%Setup learning parameters
%keep same as before
fixed_Phi_b_Trials = []; %hold all the trials
F = zeros(num_ilc_inputs, num_output_basis);

for phi_num = 1:num_ilc_inputs %try every combo possible
    %Update basis space
    phi = conjugate_basis_functions(:, phi_num);
    %Learn the controller
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q, fixed_phi_R(phi_num, phi_num), num_controllers, exploration_mag, phi, output_basis_functions, num_converged);
    F(phi_num, :) = F_phi(:, :, end); %save the controller in the stack
    fixed_Phi_b_Trials = [fixed_Phi_b_Trials, temp_Trial]; %save teh trial info
end
F_lqr = discounted_LQR(eye(num_basis_output), -pinv(output_basis_functions) * P * conjugate_basis_functions, gamma, Q, fixed_phi_R);
%Show ow they compare
sprintf('F_lqr vs Policy when Fixed Conjugate Output Basis')
lqr = F_lqr([1, 2, 19, 20], [1, 2, 19, 20])
policy = F([1, 2, 19, 20], [1, 2, 19, 20])

%Plot the path
if plot_all
    plot_dual_ilc('c Conjugate Basis', 'Scaled Conjugate Basis', fixed_Phi_b_Trials, y_star, -1, 4, update_file_path);

    F
    plot_dual_ilc('Scaled Conjugate Basis Application', 'Scaled Conjugate Basis Application', final_conj_Trial, y_star, -1, 4, update_file_path);
end

A_fixed = eye(num_ilc_states);
B_fixed = -pinv(output_basis_functions) * P * conjugate_basis_functions;
poles = eig(A_fixed + B_fixed * F); %A+BF pole formulation
mags = abs(poles);
if any(mags > 1)
    sprintf('Unstable')
    mags(mags > 1)
end

%% Growing Resolution T_y = Phi^b
growing_phi_R = 1 * eye(num_ilc_inputs);
Q = 100 * eye(num_ilc_states); %potential for full range
exploration_mag = 1; 
T_y_scale = 10;

num_converged = 0;
num_controllers = 5;
%Conjugate episodes holder
clear 'Episode';
Episode = []; 
%Setup learning parameters
%keep same as before
growing_Trials = []; %hold all the trials
F = zeros(num_ilc_inputs, num_output_basis);

for phi_num = 1:num_ilc_inputs %try every combo possible
    %Update basis space
    Episode = generate_iterative_conjugate(P, d, y_star, Episode, Q);
    phi = Episode(phi_num).Phi; %get the most recent phis
    %Learn the controller
    output_basis = phi * T_y_scale;
    input_basis = phi(:, phi_num);
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q(1:phi_num, 1:phi_num), growing_phi_R(phi_num, phi_num), num_controllers, exploration_mag, input_basis, output_basis, num_converged);
    F(phi_num, :) = [F_phi(:, :, end), zeros(1, num_ilc_inputs - phi_num)]; %save the controller in the stack
    growing_Trials = [growing_Trials, temp_Trial]; %save teh trial info
end
Phi = Episode(end).Phi;
output_basis = Phi * T_y_scale;
input_basis = Phi;
F_lqr = discounted_LQR(eye(num_basis_output), -pinv(output_basis) * P * input_basis, gamma, Q, growing_phi_R);

%Show ow they compare (bottom rows should match)
sprintf('F_lqr vs Policy when Growing Conjugate Output Basis')
lqr = F_lqr([1, 2, 19, 20], [1, 2, 19, 20]);
policy = F([1, 2, 19, 20], [1, 2, 19, 20]);

%Check the poles
A_growing = eye(num_ilc_states);
B_growing = -pinv(output_basis) * P * input_basis;
poles = eig(A_growing + B_growing * F); %A+BF pole formulation
mags = abs(poles);
if any(mags > 1)
    sprintf('Unstable')
    mags(mags > 1)
end

%Plot the path
if plot_all
    plot_dual_ilc('Growing Basis on IO', 'Growing Basis on IO', growing_Trials, y_star, -1, 4, update_file_path);

    growing_Trial = basis_ilc_sim(P, d, C*x0, input_basis, output_basis, y_star, 100, F); 
    plot_dual_ilc('Growing Basis on IO', 'Growing Basis on IO', growing_Trial, y_star, -1, 4, update_file_path);
end

%% Rolling Resolution T_y
rolling_phi_R = 1 * eye(num_ilc_inputs);
Q = 100 * eye(num_ilc_states); %potential for full range
exploration_mag = 1; 
T_y_scale = 1;

num_converged = 0;
num_controllers = 5;
%Conjugate episodes holder
clear 'Episode';
Episode = []; 
%Setup learning parameters
%keep same as before
rolling_Trials = []; %hold all the trials
F = zeros(num_ilc_inputs, num_output_basis);

for phi_num = 1:num_ilc_inputs %try every combo possible
    %Update basis space
    Episode = generate_iterative_conjugate(P, d, y_star, Episode, Q);
    phi = Episode(phi_num).Phi; %get the most recent phis
    %Learn the controller
    input_basis = phi(:, phi_num);
    output_basis = input_basis * T_y_scale;
    [temp_Trial, F_phi, ~, ~] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q(phi_num, phi_num), rolling_phi_R(phi_num, phi_num), num_controllers, exploration_mag, input_basis, output_basis, num_converged);
    F(phi_num, phi_num) = F_phi(:, :, end); %save the controller in the stack
    rolling_Trials = [rolling_Trials, temp_Trial]; %save teh trial info
end
Phi = Episode(end).Phi;
output_basis = Phi * T_y_scale;
input_basis = Phi;
F_lqr = discounted_LQR(eye(num_basis_output), -pinv(output_basis) * P * input_basis, gamma, Q, rolling_phi_R);

%Show ow they compare (bottom rows should match)
sprintf('F_lqr vs Policy when Rolling Conjugate Output Basis')
lqr = F_lqr([1, 2, 19, 20], [1, 2, 19, 20]);
policy = F([1, 2, 19, 20], [1, 2, 19, 20]);

%Check the poles
A_growing = eye(num_ilc_states);
B_growing = -pinv(output_basis) * P * input_basis;
poles = eig(A_growing + B_growing * F); %A+BF pole formulation
mags = abs(poles);
if any(mags > 1)
    sprintf('Unstable')
    mags(mags > 1)
end

%Plot the path
if plot_all
    plot_dual_ilc('Growing Basis on IO', 'Growing Basis on IO', rolling_Trials, y_star, -1, 4, update_file_path);

    growing_Trial = basis_ilc_sim(P, d, C*x0, input_basis, output_basis, y_star, 100, F); 
    plot_dual_ilc('Growing Basis on IO', 'Growing Basis on IO', growing_Trial, y_star, -1, 4, update_file_path);
end
