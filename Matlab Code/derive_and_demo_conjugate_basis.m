%Demonstration and Determination of Conjugate Basis Functions Noah Dunleavy
%Honors Thesis under direction of Professor Minh Q. Phan Thayer School of
%Engineering
clc; clear;
addpath('Saved Data', 'Functions');
setDefaultFigProp();
plot_all = true;
update_file_path = 'C:\Users\noahd\OneDrive\Desktop\Thesis\Thesis Images\Derive and Demo Conjugate Basis';  %set this to the save location if want to update figures, or -1 if not
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
p = 100;

num_ilc_states = p * num_outputs;
num_ilc_inputs = p * num_inputs;

[P, d] = P_from_ABCD(A, B, C, D, p, x0);

%% Chebyshev Generation
max_cheby = 10;
cheby_functions = generate_chebyshev(num_ilc_inputs, max_cheby);

%% Create u* and y*
beta_star = [1, 0.2, -0.3, 4, 0, 0, 0, -1, 0, 0]'; %manually set betas 
u_star = cheby_functions * beta_star;   %from basis to total space
y_star = P * u_star + d;    %from input to output

%Learning Weights
Q = 100*eye(num_ilc_states);
R = 0*eye(num_ilc_inputs);

out1_ndx = (1:2:(num_ilc_states));  %convert stacked output to individual positions
out2_ndx = (2:2:(num_ilc_states));

%% Conjugate Basis Creation (Batch)
num_conjugate_basis = max_cheby;    %how many of the cheby functions we will stimulate the system with / how many conjugate functions to make
batch_input = cheby_functions;    %for batch input (and code simplicity)
%Batch input must be full rank to ensure wronskian exists

%Generate Batch output
batch_outputs_delta = P * batch_input;  %do not include the d term, because we want the difference in outputs, which excludes d

W = batch_outputs_delta' * Q * batch_outputs_delta; %W matrix to get wronskian from

rho_batch = chol(W);   %cholesky decomposition of W to get the optimal coeffecients for the batch
T_b_batch = batch_input / (rho_batch);  %/ is same as * inv()
H_b_batch = batch_outputs_delta / (rho_batch);  %
beta_batch = H_b_batch' * Q * (y_star - d); %determined optimal weights for given basis functions (off of e_0)

%Generate final output
batch_learned_u = T_b_batch * beta_batch;
batch_learned_y = P * (batch_learned_u) + d;

%Verify conjunct condition
batch_conj_cond = T_b_batch' * (R + P' * Q * P) * T_b_batch; %should be idenitity matrix

%% Conjugate Basis Creation (Iterative)
%We already defined Q and R above We have set our input sequences through
%cheby_functions

%Initial experiments
u0 = zeros(num_ilc_inputs, 1);
y0 = d;
u1 = cheby_functions(:, 1);
y1 = P*u1 + d;

%Compute del u1 and del y1 explictly (P*del_u also works)
Episode(1).del_u = u1 - u0;
Episode(1).del_y = y1 - y0;

%Compute W
W = Episode(1).del_u' * R * Episode(1).del_u + Episode(1).del_y' * Q * Episode(1).del_y;

%rho, phi_1, h_1, and beta_1
Episode(1).rho = chol(W);
phi_1 = Episode(1).del_u * Episode(1).rho^-1; 
h_1 = Episode(1).del_y * Episode(1).rho^-1;
beta_1 = h_1' * Q * (y_star - d); %use e0 for all

Episode(1).Phi = phi_1;
Episode(1).Hb = h_1;
Episode(1).Betas = beta_1;

%Third experiment for phi_2
u2 = cheby_functions(:, 2);
y2 = P*u2 + d;

%Define all our deltas wrt u0 and y0
Episode(2).del_u = u2 - u0;
Episode(2).del_y = y2 - y0;

%Episode(2).rho, phi_2, h_2, and beta_2
Episode(2).rho(1) = 1/Episode(1).rho(1) * (Episode(1).del_u' * R * Episode(2).del_u + Episode(1).del_y' * Q * Episode(2).del_y);
Episode(2).gamma = Episode(2).rho(1)^2;
Episode(2).rho(2) = sqrt(Episode(2).del_u' * R * Episode(2).del_u + Episode(2).del_y' * Q *Episode(2).del_y - Episode(2).gamma);
phi_2 = (1/Episode(2).rho(2)) * (Episode(2).del_u - phi_1 * Episode(2).rho(1));
h_2 = 1/Episode(2).rho(2) * (Episode(2).del_y - (h_1 * Episode(2).rho(1)));
beta_2 = h_2' * Q * (y_star - d);

Episode(2).Phi = [Episode(1).Phi, phi_2];
Episode(2).Hb = [Episode(1).Hb, h_2];
Episode(2).Betas = [Episode(1).Betas; beta_2];

%Fourth experiment for phi_3
u3 = cheby_functions(:, 3);
y3 = P*u3 + d;

%Define all our deltas wrt u0 and y0
Episode(3).del_u = u3 - u0;
Episode(3).del_y = y3 - y0;

%Episode(3).rho, phi_3, h_3, and beta_3
Episode(3).rho(1) = 1/Episode(1).rho(1) * (Episode(1).del_u' * R * Episode(3).del_u + Episode(1).del_y' * Q * Episode(3).del_y);
Episode(3).rho(2) = 1/Episode(2).rho(2) * (Episode(2).del_u' * R * Episode(3).del_u + Episode(2).del_y' * Q * Episode(3).del_y - (Episode(2).rho(1)*Episode(3).rho(1)));
Episode(3).gamma = Episode(3).rho(1)^2 + Episode(3).rho(2)^2;
Episode(3).rho(3) = sqrt(Episode(3).del_u' * R * Episode(3).del_u + Episode(3).del_y' * Q *Episode(3).del_y - Episode(3).gamma);
phi_3 = (1/Episode(3).rho(3)) * (Episode(3).del_u - (phi_1 * Episode(3).rho(1) + phi_2 * Episode(3).rho(2)) );
h_3 = 1/Episode(3).rho(3) * (Episode(3).del_y - ((h_1 * Episode(3).rho(1) + h_2 * Episode(3).rho(2))));
beta_3 = h_3' * Q * (y_star - d);

Episode(3).Phi = [Episode(2).Phi, phi_3];
Episode(3).Hb = [Episode(2).Hb, h_3];
Episode(3).Betas = [Episode(2).Betas; beta_3];

for b = 3:(max_cheby-1) %for the rest of the trials 
    Episode = generate_iterative_conjugate(P, d, y_star, Episode, Q, R);
end

%% Visualize Approaches
if plot_all
    %Plot control under batch (output, we already have input%
    fig = figure('Name', 'Batched Learned Input');
    plot(batch_learned_y(out1_ndx), batch_learned_y(out2_ndx));
    hold on;
    plot(y_star(out1_ndx), y_star(out2_ndx));
    hold off;
    legend('Learned', 'Goal')
    xlabel('Mass 1 Position (m)')
    ylabel('Mass 2 Position (m)')
    title(sprintf('Batch LQL Output - 10 Conjugate Basis Functions'))
    axis equal
    save_figure(update_file_path, fig);

    %Iterative Plots Plot after 1
    learned_u = Episode(1).Phi * Episode(1).Betas;
    learned_y = P * learned_u + d;
    fig = figure('Name', 'Iterative Learned Input - 1 Basis');
    plot(learned_y(out1_ndx), learned_y(out2_ndx));
    hold on;
    plot(y_star(out1_ndx), y_star(out2_ndx));
    hold off;
    legend('Learned', 'Goal')
    xlabel('Mass 1 Position (m)')
    ylabel('Mass 2 Position (m)')
    title(('Output for LQL with 1 Conjugate Basis Function'))
    axis equal
    save_figure(update_file_path, fig);

    %Plot after 3
    learned_u = Episode(3).Phi * Episode(3).Betas;
    learned_y = P * learned_u + d;
    fig = figure('Name', 'Iterative Learned Input - 3 Basis');
    plot(learned_y(out1_ndx), learned_y(out2_ndx));
    hold on;
    plot(y_star(out1_ndx), y_star(out2_ndx));
    hold off;
    legend('Learned', 'Goal')
    xlabel('Mass 1 Position (m)')
    ylabel('Mass 2 Position (m)')
    title(('Output for LQL with 3 Conjugate Basis Functions'))
    axis equal
    save_figure(update_file_path, fig);

    %Plot after 8
    learned_u = Episode(8).Phi * Episode(8).Betas;
    learned_y = P * learned_u + d;
    fig = figure('Name', 'Iterative Learned Input - 8 Basis');
    plot(learned_y(out1_ndx), learned_y(out2_ndx));
    hold on;
    plot(y_star(out1_ndx), y_star(out2_ndx));
    hold off;
    legend('Learned', 'Goal')
    xlabel('Mass 1 Position (m)')
    ylabel('Mass 2 Position (m)')
    title(('Output for LQL with 8 Conjugate Basis Functions'))
    axis equal
    save_figure(update_file_path, fig);
end

%% Arbitrary Shape
p = 200;
loadedShape = load('Saved Data\heart_p200.mat');    %read in the file
[P_arb, d_arb] = P_from_ABCD(A, B, C, D, p, x0);
num_ilc_states = height(P_arb);
num_ilc_inputs = width(P_arb);
Q = 100 * eye(num_ilc_states);
R = 0 * eye(num_ilc_inputs);

scale = 10;
y_star_x = scale * loadedShape.drawn_x';
y_star_y = scale * loadedShape.drawn_y';
goal_matrix = [y_star_x, y_star_y]; %stack inputs next to eachother
y_star = reshape(goal_matrix', [], 1); %combine the seperate goals of each output into one vertical vector, alternating as necessary

num_to_try = 100;
[arb_conj, arb_betas] = generate_conjugate(num_ilc_inputs, num_to_try, P_arb, Q, R, d_arb, y_star);

%% Plot Arbitrary Progression
num_basis_to_include = [1, 5, 10, 20, 50, 100, num_to_try];

if plot_all
    for num_basis = num_basis_to_include
        learned_u = arb_conj(:, 1:num_basis) * arb_betas(1:num_basis);
        learned_y = P_arb * learned_u + d_arb;
        fig = figure('Name', sprintf('Iterative Learned Input for Specified Shape - %.d Basis', num_basis));
        plot(learned_y(1:2:end), learned_y(2:2:end));
        hold on;
        plot(y_star_x, y_star_y);
        hold off;
        legend('Learned', 'Goal')
        xlabel('Mass 1 Position (m)')
        ylabel('Mass 2 Position (m)')
        title(sprintf('Shaped Output for LQL with %.d Conjugate Basis Functions', num_basis))
        axis equal
        %save_figure(update_file_path, fig);
    end
end
