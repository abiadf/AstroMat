function [F,STMi] = MakeF(Z,mult)
% function [F,STMi] = MakeF(Z,mult)
% 
% This function calculates the constraint vector for multiple shooting
% problems.
%
% INPUTS:
%    Z          design variable vector, (n*n_state x 1)
%    coll       structure containing collocation parameters
%
% OUTPUTS:
%    F     constraint vector
%    STMi       3D matrix of STM between relating each subsequent pair of
%               patch points
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
BndCase = mult.BndCase;

% Convert column vector of design variables into matrix of state variables
x_ppt = reshape(Z(1:end-1),[n_state n]);

% Calculate patch point times using total time
T = Z(end);
t_ppt = linspace(0,T,n);

% Compute continuity constraints
[FC,STMi,x_fin] = Con_Cont(x_ppt,t_ppt,mult);

% Compute boundary constraints, i.e. periodicity, slack variable, etc.
switch BndCase
    case 'FixEndPtOnly'
        [FB] = Con_Bounds_FixEnd(x_ppt,x_fin,mult);
    case 'Periodicity'
        [FB] = Con_Bounds_Period(x_ppt,x_fin,mult);
end

% Assemble full constraint vector
F = [FC;FB];
