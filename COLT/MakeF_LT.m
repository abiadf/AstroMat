function [F,x0,xf,C] = MakeF_LT(Z,t_seg,t_seg_d,collmat,colt)
% function [F,x0,xf,C] = MakeF_LT(Z,t_seg,t_seg_d,collmat,colt)
% 
% This function calculates the constraint vector for collocation problems 
% involving low-thrust. 
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_seg      non-normalized times for each segment (l x (N+1)/2 x m)
%    t_seg_d    non-normalized times for each segment (l x (N-1)/2 x m)
%    collmat    structure containing constant collocation matrices
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    F     constraint vector
%    x0    state values at initial boundary node of each segment
%    xf    state values at final boundary node of each segment
%    C     matrix of polynomial coefficients (l x (N+1) x n_seg)
%
% Written by R. Pritchett, 6/07/16
% Last Update: R. Pritchett, 09/29/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_coast = colt.n_coast;

% Convert column vector of design variables to 3D matrices
[zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);

% Compute defect constraints
[FC,~,x0,xf,C] = Con_Defect(xis,uis,t_seg,t_seg_d,collmat,colt);

% Compute boundary constraints
if n_coast > 0 % if coast parameters are included
    [FB] = Con_Bounds_CSI_Coast(uis,sis,x0,xf,coast_times,colt); % NOTE: This function has not been modified to match the new Con_Bounds_VSI structure
else
    [FB] = Con_Bounds_CSI(uis,sis,x0,xf,colt);
end

% Assemble full constraint vector
F = [FC;FB];
