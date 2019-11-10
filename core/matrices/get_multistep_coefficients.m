function [coeff] = get_multistep_coefficients(step_number, method, varargin)
% Function to get the coefficients that characterize an implicit
% Adams-Moulton time marching scheme
%Input =
%          step_number: number of steps (supposed to be >=1)
%          method: update scheme used (either AM, Theta or BDF)
%          varargin: if step number is 1 and varargin is passed, it
%          represents the theta of the one-step Theta method
%Output =
%         coeff: coefficients characterizing the multistep method

assert(step_number>=1 && step_number<=6)

if step_number == 1
    if nargin > 1
        coeff = [varargin{1}, 1-varargin{1}];
    else
        coeff = [1, 0];
    end
    
elseif step_number == 2
    if strcmp(method, "AM")
        coeff = [5/12, 8/12, -1/12];
    elseif strcmp(method, "BDF")
        coeff = [2/3, -4/3, 1/3];
    else
        disp("ERROR: unrecognized multistep time marching scheme")
        msgID = 'myComponent:valueError';
        msgtext = ['Unknown value inserted for the time marching scheme. Inserted value: ', method];
        ME = MException(msgID,msgtext);
        throw(ME)
    end
        
elseif step_number == 3
    if strcmp(method, "AM")
        coeff = [9/24, 19/24, -5/24, 1/24];
    elseif strcmp(method, "BDF")
        coeff = [6/11, -18/11, 9/11, -2/11];
    else
        disp("ERROR: unrecognized multistep time marching scheme")
        msgID = 'myComponent:valueError';
        msgtext = ['Unknown value inserted for the time marching scheme. Inserted value: ', method];
        ME = MException(msgID,msgtext);
        throw(ME)
    end

elseif step_number == 4
    if strcmp(method, "AM")
        coeff = [251/720, 646/720, -264/720, 106/720, -19/720];
    elseif strcmp(method, "BDF")
        coeff = [12/25, -48/25, 36/25, -16/25, 3/25];
    else
        disp("ERROR: unrecognized multistep time marching scheme")
        msgID = 'myComponent:valueError';
        msgtext = ['Unknown value inserted for the time marching scheme. Inserted value: ', method];
        ME = MException(msgID,msgtext);
        throw(ME)
    end
   
elseif step_number == 5
    if strcmp(method, "AM")
        disp("ERROR: AM multistep method with 5 steps is not implemented")
        msgID = 'myComponent:valueError';
        msgtext = 'Invalid value inserted.';
        ME = MException(msgID,msgtext);
        throw(ME)
    elseif strcmp(method, "BDF")
        coeff = [60/137, -300/137, 300/137, -200/137, 75/137, -12/137];
    else
        disp("ERROR: unrecognized multistep time marching scheme")
        msgID = 'myComponent:valueError';
        msgtext = ['Unknown value inserted for the time marching scheme. Inserted value: ', method];
        ME = MException(msgID,msgtext);
        throw(ME)
    end
 
elseif step_number == 6
    if strcmp(method, "AM")
        disp("ERROR: AM multistep method with 6 steps is not implemented")
        msgID = 'myComponent:valueError';
        msgtext = 'Invalid value inserted.';
        ME = MException(msgID,msgtext);
        throw(ME)
    elseif strcmp(method, "BDF")
        coeff = [60/147, -360/147, 450/147, -400/147, 225/147, -72/147, 10/147];
    else
        disp("ERROR: unrecognized multistep time marching scheme")
        msgID = 'myComponent:valueError';
        msgtext = ['Unknown value inserted for the time marching scheme. Inserted value: ', method];
        ME = MException(msgID,msgtext);
        throw(ME)
    end
    
end


