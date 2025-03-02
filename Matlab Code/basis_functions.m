%Basis Function Requirements and Applications
%Noah Dunleavy
%Honors Thesis under direction of Professor Minh Q. Phan
%Thayer School of Engineering
clc; clear;
addpath('Saved Data', 'Functions');
setDefaultFigProp();
plot_all = false;
update_file_path = -1; %'C:\Users\noahd\OneDrive\Desktop\Thesis\Thesis Images\Basis Functions';  %set this to the save location if want to update figures, or -1 if not
if update_file_path ~= -1
    keyboard    %ensure we have to be very concious about ever updating all the images
end

%% Basics of Basis
max_cheby = 20; %for demo purposes, maximum number of cheby functions to generate

p = 100;   %resolution of the cheby functions
cheby_x = linspace(-1, 1, p)'; %define the cheby 'x'
cheby_functions = ones(p, max_cheby);
cheby_functions(:, 2) = cheby_x;
for ndx = 3:max_cheby       %this is the recursive cheby generation. Could be from a file, but more helpful to see
    cheby_functions(:, ndx) = 2 * cheby_x .* cheby_functions(:, ndx - 1) - cheby_functions(:, ndx - 2); %build cheby out
end

%Plot some of the core chebys
core_ndx = [1, 2, 8, 20];   %use 1 as the DC offset, then jst random
if plot_all
    core_cheby_fig = figure('Name', 'Example Cheby Functions');
    stairs(1:p, cheby_functions(:, core_ndx))
    xlabel('Chebyschev Step')
    ylabel('Amplitude')
    title('Example Chebyshev Functions')
    legend({'$T_0$', '$T_1$', '$T_7$', '$T_{19}$'}, 'Interpreter', 'latex');    %-1 off of matlab ndxs
    save_figure(update_file_path, core_cheby_fig);
end

%Plot the output signal
rng('default'); rng(10);
cheby_weights = rand_range(max_cheby, 1, -3, 3);
cheby_signal = cheby_functions * cheby_weights;
if plot_all
    fig = figure('Name', 'Example Cheby Signal');
    stairs(1:p, cheby_signal)
    xlabel('Chebyshev Step')
    ylabel('Amplitude')
    title('Example Chebyshev Composite Signal')
    save_figure(update_file_path, fig);
end
%Plot the Weights
if plot_all
    fig = figure('Name', 'Cheby Weights');
    temp_plot = plot(1:(max_cheby), cheby_weights, 'Marker', 'o', 'MarkerFaceColor', 'auto');
    xlabel('Coeffecient Number')
    ylabel('Weight')
    title('Chebyshev Coeffecient Weights')
    xlim([1, max_cheby])
    xticks(1:(max_cheby))
    save_figure(update_file_path, fig);
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
num_ilc_states = p * num_outputs;
num_ilc_inputs = p * num_inputs;

[P, d] = P_from_ABCD(A, B, C, D, p, x0);

%% Chebyshev Polynomial Creation
max_cheby = 10;
cheby_functions = generate_chebyshev(num_ilc_inputs, max_cheby);
%In theory, there could be different resolution for inputs and outputs,
%but since num_inputs = num_ouputs, we can use the same for both

%% General Parameters
%We have said p = 100
%Controller in basis_ilc_sim defaults to 0.5H+
%Built chebys out of 10 functions

%% Setup: u* in Basis Space, y* is not
T_u = cheby_functions; %cheby functions is only out of 10, so grab them all
beta_star = [1, 0.2, -0.3, 4, 0, 0, 0, -1, 0, 0]'; %manually set betas 
u_star = T_u * beta_star;   %from basis to total space
y_star = P * u_star + d;    %from input to output
T_y = T_u;  %equal - we do not want T_y to capture y*
alpha_star = pinv(T_y) * y_star;    %compute what our goal alpha is
y_Ty_alpha = T_y * alpha_star;
u_star_Tu_u_star_error = norm(y_star - y_Ty_alpha); %if fully capured, should be 0
u_star_Tu_b_star_error = norm(u_star - T_u*beta_star);  %should be zero since forced to be within

%% Render: u* in Basis Space, y* is not
if plot_all
    %u1* Plot
    fig = figure('Name', 'Goal Input 1');
    stairs(0:(p-1), u_star(1:2:end))
    title('Goal Input 1')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Partial output - Goal Input 1');
    %U2* plot
    fig = figure('Name', 'Goal Input 2');
    stairs(0:(p-1), u_star(2:2:end))
    title('Goal Input 2')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Partial output - Goal Input 2');

    %y1* Plot
    fig = figure('Name', 'Goal Output 1');
    stairs(0:(p-1), y_Ty_alpha(1:2:end))
    title('Goal Output 1 ')
    hold on;
    stairs(0:(p-1), y_star(1:2:end))
    hold off;
    legend({'$\Phi_y \alpha^\ast_1$', '$\underline{y}_1^\ast$'}, 'Interpreter', 'latex');
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Partial output - Goal Output 1');

    %y2* Plot
    fig = figure('Name', 'Goal Output 2');
    stairs(0:(p-1), y_Ty_alpha(2:2:end))
    title('Goal Output 2')
    hold on;
    stairs(0:(p-1), y_star(2:2:end))
    hold off;
    legend({'$\Phi_y \alpha^\ast_2$', '$\underline{y}_2^\ast$'}, 'Interpreter', 'latex');
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Partial output - Goal Output 2');
end

%% Application: u* in Basis Space, y* is not
num_trials = 20;
FIPO_ILC = basis_ilc_sim(P, d, C*x0, T_u, T_y, y_star, num_trials);
norm(FIPO_ILC(end).output_error) %overall error
proj_u_star = T_u * pinv(T_u' * T_u) * T_u' * u_star;   %project u* onto Tu
norm(proj_u_star - FIPO_ILC(end).input)
if plot_all
    to_plot = [1, 2, 5, 20];
    plot_ilc_coeffecients('Full Input Partial Output', 'Full Input Partial Output', FIPO_ILC, to_plot, beta_star, find(beta_star)', 5, update_file_path);
    plot_dual_ilc('Full Input Partial Output', 'Full Input Partial Output', FIPO_ILC, y_star, u_star, to_plot, update_file_path);
end 

%% Setup: y* in Basis Space, u* is not
%Keep the same Tu and Ty
alpha_star = beta_star; %new alpha_star is the old beta_star
y_star = T_y * alpha_star;
u_star = pinv(P) * (y_star - d);    %calcuate the input that gets us here
beta_star = pinv(T_u) * u_star;    %compute what our goal alpha is
u_Tu_beta = T_u * beta_star;
u_star_Tu_u_star_error = norm(u_star - u_Tu_beta) %if fully capured, should be 0

%% Render: y* in Basis Space, u* is not
if plot_all
    %u1* Plot
    fig = figure('Name', 'Goal Input 1');
    stairs(0:(p-1), u_Tu_beta(1:2:end))
    hold on;
    stairs(0:(p-1), u_star(1:2:end))
    hold off;
    legend({'$\Phi_u \beta^\ast_1$', '$\underline{u}_1^\ast$'}, 'Interpreter', 'latex'); 
    title('Goal Input 1')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Partial input Full output - Goal Input 1');
    
    %U2* plot
    fig = figure('Name', 'Goal Input 2');
    stairs(0:(p-1), u_Tu_beta(2:2:end))
    hold on;
    stairs(0:(p-1), u_star(2:2:end))
    hold off;
    legend({'$\Phi_u \beta^\ast_1$', '$\underline{u}_1^\ast$'}, 'Interpreter', 'latex'); 
    title('Goal Input 2')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Partial input Full output - Goal Input 2');
    
    %y1* Plot
    fig = figure('Name', 'Goal Output 1');
    
    stairs(0:(p-1), y_star(1:2:end))
    title('Goal Output 1')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Partial input Full output - Goal Output 1'); 
    
    %y2* Plot
    fig = figure('Name', 'Goal Output 2');
    
    stairs(0:(p-1), y_star(2:2:end))
    title('Goal Output 2')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Partial input Full output - Goal Output 2'); 
end

%% Application: y* in Basis Space, u* is not
%Same number of trials
PIFO_ILC = basis_ilc_sim(P, d, C*x0, T_u, T_y, y_star, num_trials);
norm(PIFO_ILC(end).output_error)
proj_u_star = T_u * pinv(T_u' * T_u) * T_u' * u_star;   %project u* onto Tu
norm(proj_u_star - PIFO_ILC(end).input)
if plot_all
    to_plot = [1, 2, 5, 20]; %same to_plot
    plot_ilc_coeffecients('Partial Input Full Output', 'Partial Input Full Output', PIFO_ILC, to_plot, -1, 5, find(alpha_star)', update_file_path);
    plot_dual_ilc('Partial Input Full Output', 'Partial Input Full Output', PIFO_ILC, y_star, u_star, to_plot, update_file_path);
end

%% Setup: ny < nu
%Keep the same Tu and Ty
beta_star = alpha_star; %swap back the coeffecients
u_star = T_u * beta_star;
y_star = P*u_star + d;
T_y = y_star;   %set our output basis to be pmx1

%% Render: ny < nu
if plot_all
    %Inputs already rendered

    %y1* Plot
    fig = figure('Name', 'Goal Output 1');
    
    stairs(0:(p-1), y_star(1:2:end))
    title('Goal Output 1 ')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Single output - Goal Output 1');

    %y2* Plot
    fig = figure('Name', 'Goal Output 2');
    
    stairs(0:(p-1), y_star(2:2:end))
    title('Goal Output 2')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Full input Single output - Goal Output 2');
end

%% Application: ny < nu
%Same number of trials
FISO_ILC = basis_ilc_sim(P, d, C*x0, T_u, T_y, y_star, num_trials);
norm(FISO_ILC(end).output_error)
proj_u_star = T_u * pinv(T_u' * T_u) * T_u' * u_star;   %project u* onto Tu
norm(proj_u_star - FISO_ILC(end).input)
if plot_all
    to_plot = [1, 2, 5, 20]; %same to_plot
    plot_ilc_coeffecients('Full Input Single Output', 'Full Input Single Output', FISO_ILC, to_plot, beta_star, find(beta_star)', 1, update_file_path);
    plot_dual_ilc('Full Input Single Output', 'Full Input Single Output', FISO_ILC, y_star, u_star, to_plot, update_file_path);
end

%% Setup: ny > nu
T_y = T_u; %set the basis spaces back to both be chebys
alpha_star = beta_star; %swap back the coeffecients
y_star = T_y * alpha_star;
u_star = pinv(P) * (y_star - d);
T_u = u_star;   %set our output basis to be pmx1

%% Render: ny > nu
if plot_all
    %u1* Plot
    fig = figure('Name', 'Goal Input 1');
    stairs(0:(p-1), u_star(1:2:end))
    title('Goal Input 1')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Single input Full output - Goal Input 1');
    
    %U2* plot
    fig = figure('Name', 'Goal Input 2');
    stairs(0:(p-1), u_star(2:2:end))
    title('Goal Input 2')
    %subtitle('Defined in $\Phi_u$')
    xlabel('Step Number (k)')
    ylabel('Amplitude')
    save_figure(update_file_path, fig, 'Single input Full output - Goal Input 2');

    %Outputs already rendered
end

%% Application: ny > nu
%Same number of trials
SIFO = basis_ilc_sim(P, d, C*x0, T_u, T_y, y_star, num_trials);
norm(SIFO(end).output_error)
proj_u_star = T_u * pinv(T_u' * T_u) * T_u' * u_star;   %project u* onto Tu
norm(proj_u_star - SIFO(end).input)
if plot_all
    to_plot = [1, 2, 5, 20]; %same to_plot
    plot_ilc_coeffecients('Single Input Full Output', 'Single Input Full Output', SIFO, to_plot, 1, 1, find(alpha_star)', update_file_path);
    plot_dual_ilc('Single Input Full Output', 'Single Input Full Output', SIFO, y_star, u_star, to_plot, update_file_path);
end

%% Setup: Rolling Tu
num_rolling_basis_input = 3; %num_basis - 1
num_basis_output = num_rolling_basis_input + 1; %why do anything worse than optimal
T_u_ndx = 2:(num_rolling_basis_input+1); %initilaize off a littl
T_u_full = cheby_functions; %grab max functions we will try
T_y_full = T_u_full;
T_u = T_u_full(:, [1, T_u_ndx]); %grab the first couple of basis
T_y = T_u;
%Alpha star assigned rolling
beta_star = [1, 0.2, -0.3, 4, 0, 0, 0, -1, 0, 0]'; %manually set betas 
u_star = T_u_full * beta_star;   %from basis to total space
y_star = P * u_star + d;    %from input to output

num_loops = 5;%how many loops through we'll do
num_rolls_per_loop = max_cheby / num_rolling_basis_input;
num_tries = ceil(num_loops * num_rolls_per_loop);

%% Render: Rolling Tu
if plot_all
    %Inputs already rendered

    %Outputs already rendered
end

%% Application: Rolling Tu
%Same number of trials
rolling_ILC = [];
learned_inputs = zeros(num_ilc_inputs, num_tries + 1);

for ndx = 1:num_tries
    %Run the ilc trial with basis spaces
    temporary_ILC = basis_ilc_sim(P, d, C*x0, T_u, T_y, y_star, num_trials); %sim out a trial
    rolling_ILC = [rolling_ILC, temporary_ILC]; %tack it on to previous results
    
    %What input did we learn
    learned_u = temporary_ILC(end).input;
    learned_inputs(:, ndx+1) = learned_u;
    %Take learned input, make it a bsis function
    T_u_ndx = mod(T_u_ndx + num_rolling_basis_input - 1, max_cheby) + 1; %update which new functions we'll add
    T_u = [learned_u, T_u_full(:, T_u_ndx)];
    %T_y doesnt matter
    T_y = T_u;
end

output_error = norm(rolling_ILC(end).output_error)
input_error = norm(T_u(:, 1) - u_star)

if plot_all
    to_plot = 4; 
    plot_ilc_coeffecients('Rolling Input Basis', 'Rolling Input Basis', rolling_ILC, to_plot, -1, num_rolling_basis_input + 1, num_rolling_basis_input + 1, update_file_path);
    plot_dual_ilc('Rolling Input Basis', 'Rolling Input Basis', rolling_ILC, y_star, u_star, to_plot, update_file_path);
end

