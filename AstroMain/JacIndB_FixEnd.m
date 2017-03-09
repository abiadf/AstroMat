function [iRowB,jColB] = JacIndB_FixEnd(mult)
% function [iRowB,jColB] = JacIndB_FixEnd(mult)
% 
% This function calculates the row and column indices of the nonzero 
% values of the boundary constraint Jacobian. These indices are later 
% used in the construction of the sparse DF matrix.
%
% INPUTS:
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    iRowB    row indices of the boundary constraint jacobian 
%    jColB    column indices of the boundary constraint jacobian 
%
% Originally Written by: R. Pritchett, 02/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
i_hypr0 = mult.i_hypr0; % defines which component of the initial state is constrained to a hyperplane 
x_hypr0 = mult.x_hypr0; % defines the value(s) of the initial state component(s) that define hyperplane 
n_hypr0 = mult.n_hypr0; % defines the number of initial states constrained to a hyperplane
i_hyprf = mult.i_hyprf; % defines which component of the final state is constrained to a hyperplane 
x_hyprf = mult.x_hyprf; % defines the value(s) of the final state component(s) that define hyperplane 
n_hyprf = mult.n_hyprf; % defines the number of final states constrained to a hyperplane

% Preallocate matrices
Dhypr0_iRow = zeros(n_hypr0,1);
Dhyprf_iRow = zeros(n_hyprf,1);

% Calculate number of rows taken up by continuity constraints
Cont_iRow = n_state*(n-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Matrices with Row and Column Indices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Row and column indices of hyperplane constraint for initial patch point
for ii = 1:n_hypr0
    Dhypr0_iRow(ii,1) = Cont_iRow + ii;
end
Dhypr0_jCol = i_hyprf';
    
% Row and column indices of hyperplane constraint for final patch point
for ii = 1:n_hyprf
    Dhyprf_iRow(ii,1) = Cont_iRow + n_hypr0 + ii;
end
Dhyprf_jCol = (n_state*(n-1)).*ones(n_hyprf,1) + i_hyprf';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Column Vectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRowB=[ Dhypr0_iRow;
        Dhyprf_iRow ];
        
jColB=[ Dhypr0_jCol
        Dhyprf_jCol ];
