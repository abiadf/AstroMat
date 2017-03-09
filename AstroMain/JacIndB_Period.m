function [iRowB,jColB] = JacIndB_Period(mult)
% function [iRowB,jColB] = JacIndB_Period(mult)
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
period = mult.period; % defines which states are included in periodicity constraint
n_period = mult.n_period; % defines number of states included in periodicity constraint
i_hypr0 = mult.i_hypr0; % defines which component of the initial state is constrained to a hyperplane 
n_hypr0 = mult.n_hypr0;

% Preallocate matrices
Dprdc0_iRow = zeros(n_period,1); 
Dprdcf_dx_iRow = zeros(n_period,n_state); 
Dprdcf_dt_iRow = zeros(n_period,1); 
Dhypr0_iRow = zeros(n_hypr0,1);

% Calculate number of rows taken up by continuity constraints
Cont_iRow = n_state*(n-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Matrices with Row and Column Indices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rows of periodicity constraints
for ii = 1:n_period 
    Dprdc0_iRow(ii,1) = Cont_iRow + ii;
    Dprdcf_dx_iRow(ii,:) = Cont_iRow + ii;
    Dprdcf_dt_iRow(ii,1) = Cont_iRow + ii;
end

% Columns of periodicity constraints
Dprdc0_jCol = period';
Dprdcf_dx_jCol = repmat((n_state*(n-2))+[1:n_state],[n_period 1]);
Dprdcf_dt_jCol = (n_state*n+1).*ones(n_period,1);

% Reshape Dprdcf_dx_iRow and Dprdcf_dx_jCol into column vectors
Dprdcf_dx_iRow = reshape(Dprdcf_dx_iRow,[n_period*n_state 1]);
Dprdcf_dx_jCol = reshape(Dprdcf_dx_jCol,[n_period*n_state 1]);

% Row and column indices of hyperplane constraint for initial patch point
for ii = 1:n_hypr0
    Dhypr0_iRow(ii,1) = Cont_iRow + n_period + ii;
end
Dhypr0_jCol = i_hypr0';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Column Vectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRowB=[ Dprdc0_iRow;
        Dprdcf_dx_iRow;
        Dprdcf_dt_iRow;
        Dhypr0_iRow ];
        
jColB=[ Dprdc0_jCol;
        Dprdcf_dx_jCol;
        Dprdcf_dt_jCol;
        Dhypr0_jCol ];
