function [J,gradient] = DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat)
% function [J,gradient] = DrctTrans_Obj_MaxMf(Z,t_bnd,colt)
% 
% This function calculates the value and gradient of the objective 
% funtion when the objective is to maximize the final mass of a low-thrust 
% spacecraft. The function is formatted such that it can be called directly 
% by Matlab's fmincon optimization algorithm.
%
% INPUTS:
%    Z           design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd       non-normalized times at  boundary nodes (n_seg x 1)
%    colt        structure containing collocation and optimization parameters
%    collmat    structure containing constant collocation matrices
%
% OUTPUTS:
%    J           scalar value of the objective function
%    gradient    gradient of the objective function with respect to all of 
%                the design variables
%
% Written by R. Pritchett, 6/07/16
% Last Update: R. Pritchett, 10/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Objective %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain final node states using defect constraints function
[~,~,xf] = Opt_Con_Defect(Z,colt,collmat);

% Maximize final mass value
J = -xf(7,end,end); % extract final mass value from state variable matrix
% J = 1; % extract final mass value from state variable matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Gradient of Objective Function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gradient] = DrctTrans_ObjGrad_MaxMf(Z,colt,collmat);

