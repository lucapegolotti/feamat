function [coeff] = get_multistep_coefficients(step_number, varargin)
% Function to get the coefficients that characterize an implicit
% Adams-Moulton time marching scheme
%Input =
%          step_number: number of steps (supposed to be >=1)
%          varargin: if step number is 1 and varargin is passed, it
%          represents the theta of the one-step theta method
%Output =
%         coeff: coefficients characterizing the multistep method

assert(step_number>=1 && step_number<=4)

if step_number == 1
    if nargin > 1
        coeff = [varargin{1}, 1-varargin{1}];
    else
        coeff = [1, 0];
    end
    
elseif step_number == 2
    coeff = [5/12, 8/12, -1/12];

elseif step_number == 3
    coeff = [9/24, 19/24, -5/24, 1/24];

elseif step_number == 4
    coeff = [251/720, 646/720, -264/720, 106/720, -19/720];
end


