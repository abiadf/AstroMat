function [DF_sparse_Eq,DF_sparse_InEq] = Opt_MakeDF_InEq(Z,t_bnd,ceq_l,c_l,colt)
% function [DF_sparse_InEq] = Opt_MakeDF_InEq(Z,t_bnd,colt)
% 
% This function runs setup for and then calls a function to compute the 
% partial derivatives of the constraints of the direct transcription 
% problem. It is intended to be used in functions employed by Matlab's 
% fmincon algorithm because it performs a subset of collocation related 
% operations and can be easily called by the constraint function employed 
% by fmincon.  
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd      non-normalized times at  boundary nodes (n_seg x 1)
%    ceq_l      number of equality constraints in direct transcription problem
%    c_l        number of inequality constraints in direct transcription problem
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
n_coast = colt.n_coast;
    
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
%% Sparse Matrix Setup %%
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
iRow_Eq = [iRowC; iRowB_Eq]; 
jCol_Eq = [jColC; jColB_Eq];
iRow_InEq = iRowB_InEq; 
jCol_InEq = jColB_InEq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partial Derivatives of the Jacobian %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate lengths of design and constraint vectors
Zl = length(Z);

% Compute only nonzeros of the Jacobian with Complex Step Differentiation
if n_coast > 0 % if coast parameters are included
    [DF] = MakeDF_LT_Coast(Z,t_seg,t_seg_d,collmat,colt);
else % if no coast parameters are included
    [DF_Eq] = MakeDF_LT_Eq(Z,t_seg,t_seg_d,collmat,colt);
    [DF_InEq] = MakeDF_LT_InEq(colt);
end

% Convert DF into a sparse matrix
DF_sparse_Eq = sparse(iRow_Eq,jCol_Eq,DF_Eq,ceq_l,Zl); % for complex step case
DF_sparse_InEq = sparse(iRow_InEq,jCol_InEq,DF_InEq,c_l,Zl); % for complex step case

%-------------------------------------------------------------------------%
% Finite Difference Check of Jacobian Matrix %
%-------------------------------------------------------------------------%

% % Compute full Jacobian with either Forward or Central Difference Methods
% [DF_fwrd] = MakeDF_Num_Eq(Z,ceq_l,t_seg,t_seg_d,collmat,colt);
% % DF_fwrd = zeros(ceq_l,Zl);
% 
% % Convert DF into a sparse matrix
% DF_sparse_fwrd = sparse(DF_fwrd); % for forward step case
% 
% % Convert DF_sparse into full matrix for comparison
% DF_full = full(DF_sparse_Eq);
% 
% % Calculate length of DF_sparse and DF sparsity percentage
% DFl_fwrd = nnz(DF_sparse_fwrd); % for forward step case
% DFl = length(DF_sparse_Eq); % for complex step case
% sparsity = 100*(1-DFl/(ceq_l*Zl));
% 
% % Compare forward step and sparse complex step results
% DF_comp = DF_full - DF_fwrd; % difference
% DF_comp_pcnt = (DF_comp./DF_fwrd)*100; % percent difference
% max_diff = max(DF_comp(:));
% max_pcnt_diff = max(DF_comp_pcnt(:));
% 
% % Plot forward step and complex step Jacobians (for debugging only)
% figure (1)
% hold on
% grid on
% spy(DF_sparse_Eq,'b')
% spy(DF_sparse_fwrd,'r')
% hold off

