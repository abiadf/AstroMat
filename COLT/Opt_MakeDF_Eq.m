function [DF_sparse_Eq] = Opt_MakeDF_Eq(Z,colt,collmat)
% function [DF_sparse_Eq] = Opt_MakeDF_Eq(Z,t_bnd,colt)
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
n_coast = colt.n_coast;
ceq_l = colt.ceq_l;
% c_l = colt.c_l;
    
% Extract necessary parameters from collmat stucture
t_seg = collmat.t_seg;
t_seg_d = collmat.t_seg_d;
iRow_Eq = collmat.iRow_Eq;
jCol_Eq = collmat.jCol_Eq;
% iRow_InEq = collmat.iRow_InEq;
% jCol_InEq = collmat.jCol_InEq;

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
%     [DF_InEq] = MakeDF_LT_InEq(colt);
end

% Convert DF into a sparse matrix
DF_sparse_Eq = sparse(iRow_Eq,jCol_Eq,DF_Eq,ceq_l,Zl); % for complex step case
% DF_sparse_InEq = sparse(iRow_InEq,jCol_InEq,DF_InEq,c_l,Zl); % for complex step case

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

