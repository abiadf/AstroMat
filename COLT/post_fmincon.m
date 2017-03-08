function [x_bnd,tau_nodes,C,newti,colt] = post_fmincon(Z,t_bnd,output,newti,colt)
% function [x_bnd,tau_nodes,C,newti,colt] = post_fmincon(Z,t_bnd,output,newti)
% 
% This script uses fmincon output to calculate variables needed for the
% collocation and mesh refinement process. It is intended to be run
% immediately following fmincon in the script Coll_deBoor_fmincon
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  boundary nodes (n_seg x 1)
%    output    fmincon output structure 
%    newti     Newton's method iteration counter 
%
% OUTPUTS:
%    Z          final design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    slack      column vector containing final slack variable values (n_slack x 1)
%    x_bnd      matrix of boundary node states (l x 1 x n_seg+1)
%    tau_nodes  vector of normalized variable and defect node times (N x 1) (same for all segments)
%    C          matrix of polynomial coefficients (l x (N+1) x n_seg)
%    newti      Newton's method iteration counter
%
% Written by R. Pritchett, 7/27/15
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
n_cntrl = colt.n_cntrl;
n_slack = colt.n_slack;
n_coast = colt.n_coast;
ce = colt.ce;
mu = colt.mu;
    
%-------------------------------------------------------------------------%
% Collocation Setup %
%-------------------------------------------------------------------------%

% Calculate constants
[tau_nodes,~,A_inv,B,~,~] = CollSetup(N,NodeSpace); % LG is interpolation method

% Calculate segment time intervals, divide by 2, and vectorize
t_seg = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N+1)/2])/2; % for variable nodes

% Copy constant matrices into three dimensional matrices
Ainv = repmat(A_inv,[1 1 n_seg]);
B0 = repmat(B(:,1),[1 1 n_seg]);
Bf = repmat(B(:,end),[1 1 n_seg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process fmincon Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Add fmincon Iterations to Total Iteration Count %
%-------------------------------------------------------------------------%

% Calculate number of fmincon iterations
newti = newti + output.iterations; % number of fmincon iterations; NOT Newton's Method in this case

%-------------------------------------------------------------------------%
% Break Z into 3D Matrices %
%-------------------------------------------------------------------------%

% Convert column vector of design variables to 3D matrices
[zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);

%-------------------------------------------------------------------------%
% Calculate Polynomial Coefficients %
%-------------------------------------------------------------------------%

% Calculate polynomial coefficients
fis = t_seg.*EOM_CR3BP_LTCSI_vec(xis,uis,ce,mu); % compute vector field at xis 
C = ip([xis fis],Ainv); % coefficients of polynomials

%-------------------------------------------------------------------------%
% Calculate Boundary Node Values %
%-------------------------------------------------------------------------%

% Calculate boundary node values
x0=ip(C,B0); % compute initial state of each segment
xf=ip(C,Bf); % compute final state of each segment
x_bnd = cat(3,x0,xf(:,:,end)); %Create vector of boundary nodes 

%-------------------------------------------------------------------------%
% Add fmincon output structure to colt structure %
%-------------------------------------------------------------------------%

colt.fmincon_out = output;
