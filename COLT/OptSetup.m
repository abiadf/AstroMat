function [collmat] = OptSetup(t_bnd,colt)
% function [collmat] = OptSetup(t_bnd,colt)
% 
% This function calculates variables that are constant throughout each run
% of fmincon. 
%
% INPUTS:
%    t_bnd      non-normalized times at  boundary nodes (n_seg x 1)
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    collmat    structure containing constant collocation variables
%
% Written by R. Pritchett, 03/10/17
% Last Update: R. Pritchett, 03/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
N = colt.N;
n_seg = colt.n_seg;
n_state = colt.n_state;
n_coast = colt.n_coast;
    
% Calculate constants
[~,A,A_inv,B,D,W] = CollSetup(N,NodeSpace); % LG is interpolation method

% Store constant collocation variables in single structure
collmat = struct;

% Calculate segment time intervals, divide by 2, and vectorize
collmat.t_seg = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N+1)/2])/2; % for variable nodes
collmat.t_seg_d = repmat(reshape(diff(t_bnd),[1 1 n_seg]),[n_state,(N-1)/2])/2; % for defect nodes

% Copy constant matrices into three dimensional matrices
collmat.A = repmat(A,[1 1 n_seg]);
collmat.Ainv = repmat(A_inv,[1 1 n_seg]);
collmat.Bnew = repmat(B(:,2:end-1),[1 1 n_seg]);
collmat.B0 = repmat(B(:,1),[1 1 n_seg]);
collmat.Bf = repmat(B(:,end),[1 1 n_seg]);
collmat.Dnew = repmat(D,[1 1 n_seg]);
collmat.Wnew = repmat(diag(W).',[n_state 1 n_seg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sparse Jacobian Matrix Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate indices associated with nonzeros of the collocation constraint jacobian
[iRowC,jColC] = JacIndC_LT(colt);

% Calculate indices associated with boundary constraints in the jacobian
if n_coast > 0 % if coast parameters are included
    [iRowB,jColB] = JacIndB_LT_Coast(colt);
else
    [iRowB_Eq,jColB_Eq] = JacIndB_LT_Eq(colt);
    [iRowB_InEq,jColB_InEq] = JacIndB_LT_InEq(colt);
end

% Assemble row and column indices into individual column vectors
collmat.iRow_Eq = [iRowC; iRowB_Eq]; 
collmat.jCol_Eq = [jColC; jColB_Eq];
collmat.iRow_InEq = iRowB_InEq; 
collmat.jCol_InEq = jColB_InEq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sparse Hessian Matrix Setup (Only Use if Hessian Matrix is an Input) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Calculate indices associated with nonzero collocation constraint Hessian matrices
% [jColC_hess_cell,iRowC_hess_cell] = HessIndC_LT(colt);
% 
% % Calculate indices associated with nonzero boundary constraint Hessian matrices
% [jColB_Eq_hess_cell,iRowB_Eq_hess_cell,obj_jCol_hess,obj_iRow_hess] = HessIndB_LT_Eq(colt);
% % [iRowB_InEq_hess_cell,jColB_InEq_hess_cell] = HessIndB_LT_InEq(colt);
% 
% % Combine indices associated with nonzero equality constraint Hessian matrices into single cell array
% jCol_Eq_hess_cell = cat(3,jColC_hess_cell,jColB_Eq_hess_cell);
% iRow_Eq_hess_cell = cat(3,iRowC_hess_cell,iRowB_Eq_hess_cell);
% 
% % Assemble row and column indices into individual column vectors
% collmat.jCol_Eq_hess_cell = jCol_Eq_hess_cell; 
% collmat.iRow_Eq_hess_cell = iRow_Eq_hess_cell;
% collmat.obj_jCol_hess = obj_jCol_hess; 
% collmat.obj_iRow_hess = obj_iRow_hess;
% % collmat.jColB_InEq_hess_cell = jColB_InEq_hess_cell; 
% % collmat.iRowB_InEq_hess_cell = iRowB_InEq_hess_cell; 

