function [ILC_Trial] = basis_ilc_sim(P, d, y0, T_u, T_y, y_star, num_trials, error_controller)
%Generate the ILC trials and representative data for a given system (P, d)
%with initial output y0, and basis description (T_u, T_y) to move to a goal (y_star) over specified
%trial count (num_trials)
%Inputs:
    %P: matrix - inputs u(0->(p-1)) to outputs y(1->p)
    %d: vector - noise/initial conditions matrix
    %y0: vector - initial output (since y(0) skipped by P)
    %T_u: matrix - Input basis functions (default to identity)
    %T_y: matrix - Output basis functions (default to identity)
    %y_star: vector - goal output
    %num_trials: scalar - how many trials to simulate
    %error_controller: matrix - controller of the ILC system (defaults to
    %0.8)
%Outputs:
    %ILC_Trial: structure with indexed by trial number, contains
        %inputs
        %betas (input basis weights)
        %del_betas (change in input betas)
        %outputs
        %alphas (output basis weights)
        %output error (y* - y)
        %alpha error

if T_u == -1
    T_u = eye(width(P));
end
if T_y == -1
    T_y = eye(height(P));
end

num_basis_input = width(T_u);   %basis functions are pr x n_u
num_basis_output = width(T_y);  %pm x n_y

T_y_pinv = pinv(T_y); %most common form of T_y we will use (to go from real to basis land)
alpha_star = T_y_pinv * y_star;

H = T_y_pinv * P * T_u; %descriptive IO controller, such that alpha = H * beta + T_y+ * d

if (~is_controllable(eye(num_basis_output), -H))
    fprintf('The ILC System is not Controllable with %.d basis inputs and %.d basis outputs!\n', num_basis_input, num_basis_output)
end

if ~exist('error_controller', 'var')
    error_controller = 0.5 * pinv(H);   %perfect knowledge controller as default
end

%Preallocate structure
ILC_Trial(num_trials).input = [];  %input
ILC_Trial(num_trials).betas = [];   %basis representation of input
ILC_Trial(num_trials).del_betas = [];
ILC_Trial(num_trials).output = [];  %output
ILC_Trial(num_trials).alphas = [];  %basis representation of output
ILC_Trial(num_trials).output_error = [];  %output error
ILC_Trial(num_trials).alpha_error = [];  %alpha error

%Initial trial
trial = 1;
%Inputs
ILC_Trial(trial).del_betas = zeros(num_basis_input, 1); %no 'input'
ILC_Trial(trial).betas = zeros(num_basis_input, 1);
ILC_Trial(trial).input = T_u * ILC_Trial(trial).betas;

%Outputs
relevant_output = P * ILC_Trial(trial).input + d;   %y(1) -> y(p)
ILC_Trial(trial).output = [y0; relevant_output];
ILC_Trial(trial).alphas = T_y_pinv * relevant_output; %alpha(j) = T_y+ * y(1:p)

%Errors
ILC_Trial(trial).output_error = y_star - relevant_output;
ILC_Trial(trial).alpha_error = alpha_star - ILC_Trial(trial).alphas;

trial = trial + 1;

%Simulate Rest
for trial_num = 2:num_trials
    %Inputs
    ILC_Trial(trial).del_betas = error_controller * ILC_Trial(trial - 1).alpha_error; %del_B(j) = L * e_alpha(j-1)
    ILC_Trial(trial).betas = ILC_Trial(trial - 1).betas + ILC_Trial(trial).del_betas; %B(j) = B(j-1) + del_B(j)
    ILC_Trial(trial).input = T_u * ILC_Trial(trial).betas;                            %u(j) = T_u * B(j)

    %Outputs
    relevant_output = P * ILC_Trial(trial).input + d;   %y(1) -> y(p)
    ILC_Trial(trial).output = [y0; relevant_output];
    ILC_Trial(trial).alphas = T_y_pinv * relevant_output; %alpha(j) = T_y+ * y(1:p)

    %Errors
    ILC_Trial(trial).output_error = y_star - relevant_output;
    ILC_Trial(trial).alpha_error = alpha_star - ILC_Trial(trial).alphas;

    trial = trial + 1;
end

end
