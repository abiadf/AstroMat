function [FB] = Con_Bounds_FixEnd(x_ppt,x_fin,mult)
% function [FB] = Con_Bounds_FixEnd(x_ppt,x_fin,mult)
% 
% This function calculates the boundary constraints for multiple shooting
% problems.
%
% INPUTS:
%    x_ppt      nondimensional states at patch points, (n_state x n)
%    x_fin      states at final patch point obtained from propagation
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    FB         boundary constraint vector
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
i_hypr0 = mult.i_hypr0; % defines which component of the initial state is constrained to a hyperplane 
x_hypr0 = mult.x_hypr0; % defines the value(s) of the initial state component(s) that define hyperplane 
n_hypr0 = mult.n_hypr0; % defines the number of initial states constrained to a hyperplane
i_hyprf = mult.i_hyprf; % defines which component of the final state is constrained to a hyperplane 
x_hyprf = mult.x_hyprf; % defines the value(s) of the final state component(s) that define hyperplane 
n_hyprf = mult.n_hyprf; % defines the number of final states constrained to a hyperplane

% ConB 1: Constrain initial state components to hyperplane, (may change with problem)
FB_hypr0 = x_ppt(i_hyprf,1) - x_hypr0';

% ConB 2: Constrain final state components to hyperplane, (may change with problem)
FB_hyprf = x_ppt(i_hypr0,end) - x_hyprf';

% Assemble boundary constraints into single vector
FB = [FB_hypr0;FB_hyprf];
    