function [ILC_Trial, F_policy, controller_error, F_lqr] = policy_ilc(P, d, C, D, x0, y_star, gamma, Q, R, num_controllers, exploration_mag, input_basis_functions, output_basis_functions, num_converged, existing_controller, existing_T_u, existing_T_y)
%Perform the RL policy learning on an ILC system since done so much in
%thesis
%Inputs:
    %P: matrix - inputs u(0->(p-1)) to outputs y(1->p)
    %d: vector - noise/initial conditions matrix
    %C: matrix - state to output descriptor
    %D: matrix - input to output descriptor
    %x0: vector - initial state
    %y_star: vector - goal output
    %gamma: scalar - discount factor
    %Q: matrix - cost of states (errors)
    %R: matrix - cost of inputs (change in inputs)
    %num_controllers: scalar - number of controllers to learn
    %exploration_mag: scalar - range around 0 to explore (defaults to 1)
    %input_basis_functions: matrix - basis functions on the inputs (defaults to identity)
    %output_basis_functions: matrix - basis functions on the outputs (defaults to identity) 
    %num_converged: scalar - number of trials to simulate out without exploration (defaults to 0)
%Outputs:
    %ILC_Trial: structure with indexed by trial number, contains
        %inputs
        %betas (input basis weights)
        %del_betas (change in input betas)
        %outputs
        %alphas (output basis weights)
        %output error (y* - y)
        %alpha error
    %F_policy: matrix - controller learning history
    %controller_error: scalar - normalized error from LQR
    %F_lqr: matrix - goal LQR controller

%Required Parameters to System Info
num_inputs = width(D);

%Default paramters
if ~exist('num_converged', 'var')
    num_converged = 0;  %default to no converged trials
end
if ~exist('exploration_mag', 'var')
    exploration_mag = 1;
end
if ~exist('input_basis_functions', 'var')
    input_basis_functions = eye(width(P));
end
if ~exist('output_basis_functions', 'var')
    output_basis_functions = eye(height(P));
    output_basis_functions_pinv = output_basis_functions;%save on compute time
else
    output_basis_functions_pinv = pinv(output_basis_functions);
end
num_ilc_states = width(output_basis_functions);
num_ilc_inputs = width(input_basis_functions);

%Calculate optimal controller
F_lqr = discounted_LQR(eye(num_ilc_states), -output_basis_functions_pinv * P * input_basis_functions, gamma, Q, R);

alpha_star = output_basis_functions_pinv * y_star;

%Iteration Counts
Pj_dim = num_ilc_states + num_ilc_inputs;
num_collections_per_controller = Pj_dim^2;    
total_trial_count = num_controllers * num_collections_per_controller + num_converged;

%Preallocate structure
ILC_Trial(total_trial_count).betas = [];  %input
ILC_Trial(total_trial_count).betas = [];   %basis representation of input
ILC_Trial(total_trial_count).del_beta = [];
ILC_Trial(total_trial_count).output = [];  %output
ILC_Trial(total_trial_count).alphas = [];  %basis representation of output
ILC_Trial(total_trial_count).output_error = [];  %output error
ILC_Trial(total_trial_count).alpha_error = [];  %alpha error

F_policy = zeros(num_ilc_inputs, num_ilc_states, num_controllers + 1);   %start with no controller

%Prepopulate the first trial
trial_num = 1;
ILC_Trial(trial_num).betas = zeros(num_ilc_inputs, 1);   %start with no basis guessed
ILC_Trial(trial_num).input = input_basis_functions * ILC_Trial(trial_num).betas;

ILC_Trial(trial_num).output = [C*x0; d];    %open loop response is IC and then d term
ILC_Trial(trial_num).alphas = output_basis_functions_pinv * d;
ILC_Trial(trial_num).output_error = y_star - d;    %relevant error
ILC_Trial(trial_num).alpha_error = alpha_star - ILC_Trial(trial_num).alphas;

trial_num = 2; %start at second trial now
for iteration = 1:num_controllers
    Uk_stack = zeros(num_collections_per_controller, 1);
    Xk_stack = zeros(num_collections_per_controller, (Pj_dim)^2);

    %Simulate the necessary trials
    for trial = 1:num_collections_per_controller   %number of trials to collect before updatin controller
        %ILC / Basis Controller Process
        %Beta Coeffecients
        exploration_term = rand_range(num_ilc_inputs, 1, -exploration_mag, exploration_mag);    %jiggle to learn
        ILC_Trial(trial_num).del_beta = F_policy(:, :, iteration) * ILC_Trial(trial_num - 1).alpha_error + exploration_term;
        ILC_Trial(trial_num).betas = ILC_Trial(trial_num - 1).betas + ILC_Trial(trial_num).del_beta;
        %Beta to Inputs
        ILC_Trial(trial_num).input = input_basis_functions * ILC_Trial(trial_num).betas;
        %Simulate Reality
        relevant_output = P * ILC_Trial(trial_num).input + d;   %y(1) -> y(p)
        ILC_Trial(trial_num).alphas = output_basis_functions_pinv * relevant_output;
        ILC_Trial(trial_num).output = [C*x0 + D*ILC_Trial(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
        %Calculate Error
        ILC_Trial(trial_num).output_error = y_star - relevant_output;
        ILC_Trial(trial_num).alpha_error = alpha_star - ILC_Trial(trial_num).alphas;

        %RL Translation
        state = ILC_Trial(trial_num - 1).alpha_error; %analogous state
        input = ILC_Trial(trial_num).del_beta;        %analogous input
        next_state = ILC_Trial(trial_num).alpha_error;    %x(k+1) = e_alpha(j)
        next_input = F_policy(:, :, iteration) * next_state; %no exploration term here

        xu_stack = [state; input];
        xu_next_stack = [next_state; next_input];

        Xk_stack(trial, :) = kron(xu_stack', xu_stack') - gamma * kron(xu_next_stack', xu_next_stack');
        Uk_stack(trial, :) = input' * R * input + state' * Q * state;

        trial_num = trial_num + 1;
    end

    %Calculate P and new controller
    PjS = pinv(Xk_stack) * Uk_stack;
    Pj = reshape(PjS, Pj_dim, Pj_dim);
    Pj = 0.5 * (Pj + Pj'); %to impose symmetry (significantly reduces error)
    Pjuu = Pj((num_ilc_states+1):end, (num_ilc_states+1):end);
    PjxuT = Pj((num_ilc_states+1):end, 1:num_ilc_states);
    new_F = -pinv(Pjuu) * PjxuT;
    F_policy(:, :, iteration + 1) = new_F;
end

controller_error = norm(F_policy(:, :, end) - F_lqr)/numel(F_lqr);


for ndx = 1:num_converged
    %ILC / Basis Controller Process
    %Beta Coeffecients
    ILC_Trial(trial_num).del_beta = F_policy(:, :, iteration) * ILC_Trial(trial_num - 1).alpha_error;
    ILC_Trial(trial_num).betas = ILC_Trial(trial_num - 1).betas + ILC_Trial(trial_num).del_beta;
    %Beta to Inputs
    ILC_Trial(trial_num).input = input_basis_functions * ILC_Trial(trial_num).betas;
    %Simulate Reality
    relevant_output = P * ILC_Trial(trial_num).input + d;   %y(1) -> y(p)
    ILC_Trial(trial_num).alphas = output_basis_functions_pinv * relevant_output;
    ILC_Trial(trial_num).output = [C*x0 + D*ILC_Trial(trial_num).input(1:num_inputs); relevant_output]; %total output y(0) -> y(p) for completeness
    %Calculate Error
    ILC_Trial(trial_num).output_error = y_star - relevant_output;
    ILC_Trial(trial_num).alpha_error = alpha_star - ILC_Trial(trial_num).alphas;

    trial_num = trial_num + 1;
end

end