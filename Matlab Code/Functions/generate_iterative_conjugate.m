function [Episode] = generate_iterative_conjugate(P, d, y_star, old_episodes, Q, R)
%From an exisiting conjugate basis space, apply a new episode to generate a
%new basis
%Input:
    %P: matrix -descriptie matrix to map u to y (y = Pu + d)
    %d: vector - handle noise in the relation to map u to y
    %y_star: vector - goal output
    %old_episodes: struct - structure holding old iterative data
    %Q: matrix - cost of states
    %R: matrix - cost of inputs
%Output:
    %Episode: struct - all the previous trials and info, plus the new one

num_ilc_inputs = width(P);
num_ilc_states = height(P);

if nargin < 6%if no R
    R = zeros(num_ilc_inputs);
end
if nargin < 5%if no Q
    Q = 100 * eye(num_ilc_states);
end
if ((nargin < 4) || isempty(old_episodes)) %if we have no episodes to start, do the first trial
    Episode(1).del_u = ones(num_ilc_inputs, 1); %first chebyshev is all ones
    Episode(1).del_y = P * Episode(1).del_u;

    %Compute W
    W = Episode(1).del_u' * R * Episode(1).del_u + Episode(1).del_y' * Q * Episode(1).del_y;

    %rho_1, phi_1, h_1, and beta_1
    Episode(1).rho = chol(W);
    phi_1 = Episode(1).del_u * Episode(1).rho^-1;
    h_1 = Episode(1).del_y * Episode(1).rho^-1;
    beta_1 = h_1' * Q * (y_star - d); %use e0 for all

    Episode(1).Phi = phi_1;
    Episode(1).Hb = h_1;
    Episode(1).Betas = beta_1;
    return
end

%If this is not our first trial
Episode = old_episodes;

b = length(Episode);

%Generate our dels
tried_inputs = zeros(num_ilc_inputs, b);
for ndx = 1:b
    tried_inputs(:, ndx) = Episode(ndx).del_u;
end

gen_cheby = generate_chebyshev(num_ilc_inputs, b + 1); %there is likely a more effecient way, but to ensure all the generations are exactly the same
next_input = gen_cheby(:, b+1); %next input is our final cheby
while (rank([tried_inputs, next_input]) < (b+1))    %ensure that our added input does not ruin the rank of the system
    next_input = rand_range(num_ilc_inputs, 1, -1, 1);
    return
end

Episode(b + 1).del_u = next_input;   %-u0, but open loop
Episode(b + 1).del_y = P * Episode(b + 1).del_u;

%Compute first component
Episode(b + 1).rho(1) = 1/Episode(1).rho(1) * (Episode(1).del_u' * R * Episode(b + 1).del_u + Episode(1).del_y' * Q * Episode(b + 1).del_y);
%Middle Components
if b >= 2
    for i = 2:b
        sum_term = 0;
        for j = 1:(i-1)
            sum_term = sum_term + Episode(i).rho(j) * Episode(b + 1).rho(j); %hard to vector format across structs
        end
        Episode(b+1).rho(i) = (1 / Episode(i).rho(i)) * (Episode(i).del_u' * R * Episode(b+1).del_u + Episode(i).del_y' * Q * Episode(b+1).del_y - sum_term);
    end
end

%Last Component
Episode(b + 1).gamma = Episode(b + 1).rho(1:b) * Episode(b + 1).rho(1:b)';
Episode(b + 1).rho(b+1) = sqrt(Episode(b + 1).del_u' * R * Episode(b + 1).del_u + Episode(b + 1).del_y' * Q *Episode(b + 1).del_y - Episode(b + 1).gamma);

%Phi calculation
new_basis = (1/Episode(b + 1).rho(b+1)) * (Episode(b + 1).del_u - Episode(b).Phi * Episode(b + 1).rho(1:b)');
Episode(b+1).Phi = [Episode(b).Phi, new_basis];
%Hb
new_h = (1/Episode(b + 1).rho(b+1)) * (Episode(b + 1).del_y - Episode(b).Hb * Episode(b + 1).rho(1:b)');
Episode(b+1).Hb = [Episode(b).Hb, new_h];
%Betas
new_beta = new_h' * Q * (y_star - d);
Episode(b+1).Betas = [Episode(b).Betas; new_beta];

end


