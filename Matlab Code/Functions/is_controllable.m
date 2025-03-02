function [controllable] = is_controllable(A, B, precision)
%Check whether or not a matrix specificed by A and B matricies is
%controllable
%Inputs:
    %A: matrix - state dynamics
    %B: matrix - input dynamics
    %precision: scalar - set how precise matlab is with its rounding
%Outputs:
    %controllable: bool - whether or not the system is controllable

if exist('precision', 'var') %ability to set precision
    A = vpa(A, precision);  
    B = vpa(B, precision);
end

controllability_matrix = ctrb(A, B);
controllable = (rank(controllability_matrix) == height(A));

end