function [c,ceq,Dc,Dceq] = DrctTrans_Con(Z,t_bnd,colt,collmat)
% function [c,ceq,Dc,Dceq] = DrctTrans_Con(Z,t_bnd,colt)
% 
% This function calculates the equality and inequality constraints of the 
% low-thrust optimization problem that employs direct transcription. The 
% function is formatted such that it can be called directly by Matlab's 
% fmincon optimization algorithm.
%
% INPUTS:
%    Z           design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd       non-normalized times at  boundary nodes (n_seg x 1)
%    colt        structure containing collocation and optimization parameters
%
% OUTPUTS:
%    c           vector of inequality constraints
%    ceq         vector of equality constraints
%    Dc          matrix of partial derivatives of the inequality 
%                constraints with respect to the design variables 
%    Dceq        matrix of partial derivatives of the equality 
%                constraints with respect to the design variables
%
% Written by R. Pritchett, 03/22/2017
% Last Update: R. Pritchett, 03/22/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_coast = colt.n_coast;

% Convert column vector of design variables to 3D matrices
[~,~,uis,~,coast_times,~] = Z23D(Z,colt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Equality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute defect constraints
[FC,x0,xf] = Opt_Con_Defect(Z,t_bnd,colt,collmat);

% Compute boundary equality constraints
if n_coast > 0 % if coast parameters are included
    [FB] = Con_Bounds_VSI_Coast(uis,sis,x0,xf,coast_times,colt);
else
    [FB_eq] = Con_Bounds_CSI_Eq(uis,x0,xf,colt);
end

% Assemble equality constraints
ceq = [FC;FB_eq]; % equality constraints
colt.ceq_l = length(ceq); % number of equality constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Inequality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Compute boundary inequality constraints
% [FB_ineq] = Con_Bounds_CSI_InEq(uis,colt);
% 
% % Assemble inequality constraints
% c = FB_ineq;

% Define c as empty because this function does not inlude inequality constraints
c = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Gradient of Equality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute partial derivatives with respect to equality constraints
[Dceq] = Opt_MakeDF_Eq(Z,t_bnd,colt,collmat);
Dceq = transpose(Dceq); % Matlab requires the transpose of what my algorithm produces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Gradient of Inequality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dc = transpose(Dc); % Matlab requires the transpose of what my algorithm produces

% Define Dc as empty because this function does not inlude inequality constraints
Dc = [];