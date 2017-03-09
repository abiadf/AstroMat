function [iRowC,jColC] = JacIndC(mult)
% function [iRowC,jColC] = JacIndC(mult)
% 
% This function calculates the row and column indices of the nonzero 
% values of the continuity constraint Jacobian. These indices are later 
% used in the construction of the sparse DF matrix. Note that the output 
% column vectors do not include indices for the boundary constraints.
% Indices for the boundary constraints are calculated in JacIndB.
%
% INPUTS:
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    iRowC    row indices of the continuity constraint jacobian 
%    jColC    column indices of the continuity constraint jacobian 
%
% Originally Written by: R. Pritchett, 02/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;

% Preallocate matrices
DSTM_iRow = zeros(n_state,n_state,n-1); 
DSTM_jCol = zeros(n_state,n_state,n-1);
Dneye_iRow = zeros(n_state,n_state,n-1);
Dneye_jCol = zeros(n_state,n_state,n-1);
Ddt_iRow = zeros(n_state,n-1);
Ddt_jCol = zeros(n_state,n-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Matrices with Row and Column Indices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:n-1
    
    for jj = 1:n_state
        
        DSTM_iRow(jj,:,ii) = n_state*(ii-1) + jj;
        Dneye_iRow(jj,:,ii) = n_state*(ii-1) + jj;
        Ddt_iRow(jj,ii) = n_state*(ii-1) + jj;
        
        DSTM_jCol(:,jj,ii) = n_state*(ii-1) + jj;
        Dneye_jCol(:,jj,ii) = n_state*ii + jj;
        Ddt_jCol(jj,ii) = n_state*n + 1;
        
    end
    
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Column Vectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRowC=[ reshape(DSTM_iRow,[n_state*n_state*(n-1) 1]);
        reshape(Dneye_iRow,[n_state*n_state*(n-1) 1]);
        reshape(Ddt_iRow,[n_state*(n-1) 1]) ];
        
jColC=[ reshape(DSTM_jCol,[n_state*n_state*(n-1) 1]);
        reshape(Dneye_jCol,[n_state*n_state*(n-1) 1]);
        reshape(Ddt_jCol,[n_state*(n-1) 1]) ];
    