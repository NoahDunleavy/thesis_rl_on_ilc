function [conjugate_basis_functions, conjugate_betas] = generate_conjugate(basis_resolution,num_basis, P, Q, R, d, y_star)
%Create conjugate basis functions for a system, derived from chebyshevs
%Input:
    %height: scalar - height/resolution of the conjugate functions to generate
    %num_basis: scalar - num_basis/number of functions to generate
    %P: matrix - system matrix relating 
    %Q: matrix - cost of states
    %d: vector - IC vector
    %y_star: vector - goal output
%Output:
    %conjugate_basis_functions: matrix - funcions that satify conjunctionality for the system P
    %conjuagte_betas: vector - optimal weighting to minimize error

if nargin < 7
    y_star = zeros(height(P), 1); %default goal to zeros
end
if nargin < 6
    d = zeros(height(P), 1); %default noise to zeros
end
if nargin < 5
    R = 0 * eye(height(P)); %default input cost to 0
end
if nargin < 4
    Q = 100 * eye(width(P)); %default state costs to 100
end

%Create chebys
batch_input = generate_chebyshev(basis_resolution, num_basis);

%Generate Batch output
batch_outputs_delta = P * batch_input;  %do not include the d term, because we want the difference in outputs, which excludes d

W = batch_input' * R * batch_input + batch_outputs_delta' * Q * batch_outputs_delta; %W matrix 

while rank(W) < min(basis_resolution, num_basis) %if dont get a full rank W
    batch_input = rand_range(basis_resolution, num_basis, -10, 10); %go different input approach
    
    %Generate Batch output
    batch_outputs_delta = P * batch_input;  %do not include the d term, because we want the difference in outputs, which excludes d
    
    W = batch_input' * R * batch_input + batch_outputs_delta' * Q * batch_outputs_delta; %W matrix 
    sprintf('Using random inputs')
end

rho_batch = chol(W);   %cholesky decomposition of W to get the optimal coeffecients for the batch
T_b = batch_input / (rho_batch);  %/ is same as * inv()

conjugate_basis_functions = T_b;

H_b = batch_outputs_delta / (rho_batch);
conjugate_betas = H_b' * Q * (y_star - d); %determined optimal weights for given basis functions (off of e_0)


end