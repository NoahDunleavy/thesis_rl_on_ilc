function [num] = rand_range(height, width, lower, upper)
%Return a random number in a range
%Inputs:
    %height: scalar - height (num rows) of random vector
    %width: scalar - width (num columns)
    %lower: scalar - lower bound of numbers to generate
    %upper: scalar - upper bound
%Outputs:
    %num: matrix - height x width matrix of random numbers in range

num = lower + rand(height, width)*(upper - lower); %random times the amplitude, shift by lower bound

end