function [FB] = Con_Bounds_Period(x_ppt,x_fin,mult)
% function [FB] = Con_Bounds_Period(x_ppt,mult)
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
period = mult.period; % defines which states are included in periodicity constraint
i_hypr0 = mult.i_hypr0; % defines which component of the initial state is constrained to a hyperplane 
x_hypr0 = mult.x_hypr0; % defines the value(s) of the component(s) that define the hyperplane 

% ConB 1: Compute periodicity constraint
FB_prdc = x_fin(1,period)' - x_ppt(period,1);

% ConB 2: Constrain initial y component to hyperplane, (may change with problem)
FB_hypr0 = x_ppt(i_hypr0,1) - x_hypr0';

% Assemble boundary constraints into single vector
FB = [FB_prdc;FB_hypr0];
    