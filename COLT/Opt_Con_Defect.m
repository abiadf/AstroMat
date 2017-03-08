function [FC,x0,xf] = Opt_Con_Defect(Z,t_bnd,colt)
% function [FC,x0,xf] = Opt_Con_Defect(Z,t_bnd,colt)
% 
% This function runs setup for and then calls a function to compute the 
% defect and continuity constraints associated with collocation. It is 
% intended to be used in functions employed by Matlab's fmincon algorithm 
% because it performs a subset of collocation related operations and can be 
% easily called by the objective or constraint function employed by fmincon.  
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd      non-normalized times at  boundary nodes (n_seg x 1)
%    colt       structure containing collocation and optimization parameters
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
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
N = colt.N;
n_seg = colt.n_seg;
n_state = colt.n_state;
    
% Calculate constants
[~,A,A_inv,B,D,W] = CollSetup(N,NodeSpace); % LG is interpolation method

% Calculate segment time intervals, divide by 2, and vectorize
t_seg = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N+1)/2])/2; % for variable nodes
t_seg_d = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N-1)/2])/2; % for defect nodes

% Store collocation matrices in structure
collmat = struct;

% Copy constant matrices into three dimensional matrices
collmat.A = repmat(A,[1 1 n_seg]);
collmat.Ainv = repmat(A_inv,[1 1 n_seg]);
collmat.Bnew = repmat(B(:,2:end-1),[1 1 n_seg]);
collmat.B0 = repmat(B(:,1),[1 1 n_seg]);
collmat.Bf = repmat(B(:,end),[1 1 n_seg]);
collmat.Dnew = repmat(D,[1 1 n_seg]);
collmat.Wnew = repmat(diag(W).',[n_state 1 n_seg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Defect Constraints and Boundary Nodes %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert column vector of design variables to 3D matrices
[~,xis,uis,~,~,~] = Z23D(Z,colt);

% Compute defect constraints
[FC,~,x0,xf,~] = Con_Defect(xis,uis,t_seg,t_seg_d,collmat,colt);


