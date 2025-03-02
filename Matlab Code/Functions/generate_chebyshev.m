function [cheby_functions] = generate_chebyshev(cheb_resolution, num_cheby)
%Generate a matrix of 'num_cheby' functions, with depth/resolution of
%'cheb_resolution'
%Input:
    %cheb_resolution: scalar - resolution of each functions, 'height' of matrix
    %num_cheby: scalar - number of chebyshevs to generate
%Output:
    %cheby_functions: matrix - chebyshev functions, type 1

cheby_x = linspace(-1, 1, cheb_resolution)'; %define the cheby 'x'

cheby_functions = ones(cheb_resolution, num_cheby);
cheby_functions(:, 2) = cheby_x;

for ndx = 3:num_cheby       %this is the recursive cheby generation. Could be from a file, but more helpful to see
    cheby_functions(:, ndx) = 2 * cheby_x .* cheby_functions(:, ndx - 1) - cheby_functions(:, ndx - 2); %build cheby out
end

end