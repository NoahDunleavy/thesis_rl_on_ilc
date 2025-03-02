function [F] = discounted_LQR(A, B, discount_factor, Q, R, verbose)
%From system specifications, determine the optimal feedback controller
%(Linear quadratic regulator). Of form u(k) = F*x(k) optimally
%Inputs:
    %A: matrix - system dynamics
    %B: matrix - input dynamics
    %discount_factor: scalar - discount factor
    %Q: matrix - state costs
    %R: matrix - input costs
    %verbose: bool - system dynamics
%Outputs:
    %F: matrix - found controller

%Check if the verbose provided
if nargin < 6
    verbose = false; %Default false
end

%Solve for P - the solution to Algebraic Riccati Equation (Discrete - DARE):
R_gamma = R/discount_factor;
A_gamma = sqrt(discount_factor)*A;
[P, ~, ~] = idare(A_gamma, B, Q, R_gamma); %solve the riccati equation

%Verify it is a solution (text output for user)
if verbose
    test_P = A_gamma'*P*A_gamma - A_gamma'*P*B / ((R_gamma + B'*P*B)) * B'*P*A_gamma + Q; %verify P is a solution
    sprintf('Solved Ps suitability to solve the Riccati equation is %g', norm(test_P - P)) %names the answer this way - want this small
end

%Calculate F
F = (-1/sqrt(discount_factor)) * inv((transpose(B) * P * B + R_gamma)) * transpose(B)*P*A_gamma; %LQR solution
end