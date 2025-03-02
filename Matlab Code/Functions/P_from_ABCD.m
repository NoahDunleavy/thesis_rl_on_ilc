function [P, d] = P_from_ABCD(A, B, C, D, p, x0)
%Construct the P matrix and d matrix for given matrix values
%Will satisfy the equations y_bar = P * u_bar + d
%y is y(1) -> y(p)
%u is u(0) -> u(p-1)
%d captures initial conditions and noise
%Inputs:
    %A: matrix - state dynamics matrix
    %B: matrix - input dynamics matrix
    %C: matrix - state to output
    %D: matrix - input to output
    %p: scalar - number of steps to map out
    %x0: vector - initial conditions
%Outputs:
    %P: matrix - ILC system mapping matrix for u(0->(p-1)) -> y(1->p)
    %d: vector - captures the 'disturbance' caused by initial conditions

num_states = width(A);
num_inputs = width(B);
num_outputs = height(C);

%Construct 'True' P
P = zeros(num_states, num_inputs);
ic_matrix = zeros(num_outputs, num_states);  %matrix which governs the impact of the Ics (for sim)
for row = 1:p
    row_start = ((row - 1) * num_outputs) + 1;
    row_end = row_start + num_outputs - 1;
    for col = 1:row
        col_start = ((col - 1) * num_inputs) + 1;
        col_end = col_start + num_inputs - 1;
        if (row + 1 == col)
            P(row_start:row_end, col_start:col_end) = D;
        else
            mat_pow = row - col;
            P(row_start:row_end, col_start:col_end) = C * (A^mat_pow) * B;
        end

    end
    ic_matrix(row_start:row_end, :) = C * A^row;  %construct the IC matrix
end
d = ic_matrix * x0;

%Verify P is accurate
demo_in = rand_range(num_inputs * p, 1, -2, 2);
demo_P_out = P * demo_in + d;   %y = Pu + d

demo_dlsim_in = reshape(demo_in, num_inputs, [])';  %dlsim takes a pxr matrix, whereas P is a pr x 1
demo_dlsim_out_matrix = dlsim(A, B, C, D, [demo_dlsim_in; zeros(1, num_inputs)], x0);
demo_dlsim_out = reshape(demo_dlsim_out_matrix(2:end, :)', [], 1);

if (norm(demo_P_out - demo_dlsim_out)/numel(demo_P_out) > 1e-6)
    fprintf('P does not capture output accurately')
end


end