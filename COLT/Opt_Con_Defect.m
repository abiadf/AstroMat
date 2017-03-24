function [FC,x0,xf] = Opt_Con_Defect(Z,colt,collmat)
% function [FC,x0,xf] = Opt_Con_Defect(Z,colt,collmat)
% 
% This function runs setup for and then calls a function to compute the 
% defect and continuity constraints associated with collocation. It is 
% intended to be used in functions employed by Matlab's fmincon algorithm 
% because it performs a subset of collocation related operations and can be 
% easily called by the objective or constraint function employed by fmincon.  
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    colt       structure containing collocation and optimization parameters
%    collmat    structure containing constant collocation matrices
%
% OUTPUTS:
%    FC     constraint vector, omitting boundary constraints
%    x0     initial state for segment  (l x 1 x n_seg)
%    xf     final state for segment  (l x 1 x n_seg)
%
% Written by R. Pritchett, 6/07/16
% Last Update: R. Pritchett, 10/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Defect Constraints and Boundary Nodes %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from collmat stucture
t_seg = collmat.t_seg;
t_seg_d = collmat.t_seg_d;

% Convert column vector of design variables to 3D matrices
[~,xis,uis,~,~,~] = Z23D(Z,colt);

% Compute defect constraints
[FC,~,x0,xf,~] = Con_Defect(xis,uis,t_seg,t_seg_d,collmat,colt);

