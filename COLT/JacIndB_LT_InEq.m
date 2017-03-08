function [iRowB_InEq,jColB_InEq] = JacIndB_LT_InEq(colt)
% function [iRowB_InEq,jColB_InEq] = JacIndB_LT_InEq(colt)
% 
% This function calculates the row and column indices of the nonzero 
% values of the boundary inequality constraints in the collocation constraint 
% jacobian when coast parameters are NOT included. These indices are later used
% in the construction of the sparse DF matrix.
%
% INPUTS:
%    colt     structure containing collocation and optimization parameters
%    
% OUTPUTS:
%    iRowB_InEq   row indices of the boundary constraints in the collocation inequality constraint jacobian 
%    jColB_InEq   column indices of the boundary constraints in the collocation inequality constraint jacobian  
%
% Written by R. Pritchett, 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N; % degree of interpolating polynomial
n_seg = colt.n_seg; % number of segments
n_state = colt.n_state; % number of states, same as number of equations of motion
n_cntrl = colt.n_cntrl; % number of control variables
n_slack = colt.n_slack; % number of slack variables

%% Preallocate Matrices %%

% Total design variables per variable node
l = (N+1)*n_state/2 + n_cntrl + n_slack;

%Preallocate matrices of partial derivatives for boundary and slack constraints
D_ineq_iRow = zeros(2,n_seg);
D_ineq_jCol = zeros(2,n_seg);
    
%% Fill Matrices of InEquality Constraints %%

% Fill row number matrices
D_ineq_iRow(1,:) = [1:2:2*n_seg];
D_ineq_iRow(2,:) = [2:2:2*n_seg];

% Fill column number matrices
for ii = 1:n_seg
    D_ineq_jCol(:,ii) = (ii-1)*l + 1;
end

%% Assemble Column Vectors %%
    
iRowB_InEq = reshape(D_ineq_iRow,[2*n_seg 1]);
      
jColB_InEq = reshape(D_ineq_jCol,[2*n_seg 1]);
    