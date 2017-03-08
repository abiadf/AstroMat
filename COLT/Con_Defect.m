function [FC,delta,x0,xf,C] = Con_Defect(xis,uis,t_seg,t_seg_d,collmat,colt)
% function [FC,delta,x0,xf,C] = Con_Defect(xis,uis,t_seg,t_seg_d,collmat,colt)
% 
% This function computes the constraint vector, omitting boundary 
% constraints, for the collocation problem. This constraint vector includes
% continuity constraints and defect node constraints.
%
% INPUTS:
%    xis        design variable matrix, (n_state*n_seg*(N+1)/2 x 1)
%    uis        control variable matrix (n_cntrl*n_seg*(N+1)/2 x 1)
%    t_seg      non-normalized times for each segment (l x (N+1)/2 x m)
%    t_seg_d    non-normalized times for each segment (l x (N-1)/2 x m)
%    collmat    structure containing constant collocation matrices
%    colt       structure containing collocation and optimization parameters
% 
% OUTPUTS:
%    FC     constraint vector, omitting boundary constraints
%    delta  collocation defect constraints (l x (N-1)/2 x n_seg)
%    x0     initial state for segment  (l x 1 x n_seg)
%    xf     final state for segment  (l x 1 x n_seg)
%    C      matrix of polynomial coefficients
%
% See Grebow and Pavlak (2015) for details on the algorithm and more
% information about the matricies Ainv, B, and D.
%
% Originally Written by: R. Pritchett, 07/15/2015
% Last Update: R. Pritchett, 09/29/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N; % degree of collocation polynomial
n_seg = colt.n_seg; % number of segments
n_state = colt.n_state; % number of states, same as number of equations of motion
ce = colt.ce; % s/c effective exhaust velocity
mu = colt.mu; % CR3BP mass ratio

% Extract necessary parameters from collocation matrix structure
Ainv = collmat.Ainv; % inverse of full time matrix for nodes that are variables (n+1 x n+1 x m)
B = collmat.Bnew; % time matrix for nodes that are constraints (assume does not contain endpoints of segment)  (n+1 x (n-1)/2 x m)
B0 = collmat.B0; % time matrix for initial state of segment (N+1 x 1 x n_seg)
Bf = collmat.Bf; % time matrix for final state of segment (N+1 x 1 x n_seg)
D = collmat.Dnew; % derivative time matrix at constraints points (n+1 x (n-1)/2 x m)
W = collmat.Wnew; % quadrature weights associated with LG points (l x (n-1/2) x m)

% compute vector field at xis 
fis=t_seg.*EOM_CR3BP_LTCSI_vec(xis,uis,ce,mu);

% coefficients of polynomials
C=ip([xis fis],Ainv);

% interpolate to obtain states at defect nodes
xis_d=ip(C,B);

%%% Note: control variables are assumed constant across segment %%%
uis_d = uis(:,1:(N-1)/2,:); % define control variables at defect nodes

% compute vector field at defect states and times
fis_m=t_seg_d.*EOM_CR3BP_LTCSI_vec(xis_d,uis_d,ce,mu);

% defect constraints
delta=(ip(C,D)-fis_m).*W;

% compute initial and final state of each segment
x0=ip(C,B0);
xf=ip(C,Bf);

% calculate constraints and assemble into a single column vector
FC = [reshape(x0(:,:,2:end)-xf(:,:,1:end-1),(n_seg-1)*n_state,1);
      reshape(delta,n_seg*n_state*(N-1)/2,1)];
