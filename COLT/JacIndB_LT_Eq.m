function [iRowB_Eq,jColB_Eq] = JacIndB_LT_Eq(colt)
% function [iRowB_Eq,jColB_Eq] = JacIndB_LT_Eq(colt)
% 
% This function calculates the row and column indices of the nonzero 
% values of the boundary equality constraints in the collocation constraint 
% jacobian when coast parameters are NOT included. These indices are later 
% used in the construction of the sparse DF matrix for equality
% constraints.
%
% INPUTS:
%    colt     structure containing collocation and optimization parameters
%    
% OUTPUTS:
%    iRowB   row indices of the boundary constraints in the collocation constraint jacobian 
%    jColB   column indices of the boundary constraints in the collocation constraint jacobian  
%
% Written by R. Pritchett, 02/01/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N; % degree of interpolating polynomial
n_seg = colt.n_seg; % number of segments
n_state = colt.n_state; % number of states, same as number of equations of motion
n_cntrl = colt.n_cntrl; % number of control variables
n_slack = colt.n_slack; % number of slack variables
x0_des_ind = colt.x0_des_ind;
xf_des_ind = colt.xf_des_ind;

% Calculate number of fixed components in initial and final endpoints
n_fix0 = length(x0_des_ind);
n_fixf = length(xf_des_ind);

%% Preallocate Matrices %%

% Total design variables per variable node
l = (N+1)*n_state/2 + n_cntrl + n_slack;

%Preallocate matrices of partial derivatives for boundary and slack constraints
Dbnd0_iRow = zeros(n_fix0,l); % derivatives of initial endpoint constraint
Dbnd0_jCol = zeros(n_fix0,l); % derivatives of initial endpoint constraint
Dbndf_iRow = zeros(n_fixf,l); % derivatives of final endpoint constraint
Dbndf_jCol = zeros(n_fixf,l); % derivatives of final endpoint constraint
Duhat_iRow = zeros(1,3,n_seg); % derivatives of uhat magnitude constraint
Duhat_jCol = zeros(1,3,n_seg); % derivatives of uhat magnitude constraint

% Calculate number of rows taken up by collocation constraints
CollCon_iRow = n_state*(n_seg-1) + n_state*(N-1)/2*n_seg;

%% Fill Matrices of Fixed Endpoint Constraints %%

% Fill rows of initial endpoint constraint partials
for ii = 1:n_fix0
    Dbnd0_iRow(ii,:) = CollCon_iRow + ii;
end

% Fill rows of final endpoint constraint partials
for ii = 1:n_fixf
    Dbndf_iRow(ii,:) = CollCon_iRow + n_fix0 + ii;
end

% Fill columns of initial and final endpoint constraint partials
for jj = 1:l
    Dbnd0_jCol(:,jj) = jj; % initial endpoint constraint
    Dbndf_jCol(:,jj) = l*(n_seg-1) + jj; % final endpoint constraint
end
    
%% Fill Matrices of uhat Magnitude Constraints %%

%Boundary counstraints with respect to state variables
for ii = 1:n_seg 
    
    % Iterate rows of uhat magnitude constraint
    for jj = 1:3
        Duhat_iRow(1,jj,ii) = CollCon_iRow + n_fix0 + n_fixf + ii;
    end
    
    % Iterate columns of uhat magnitude constraint
    for jj = 1:3
        Duhat_jCol(1,jj,ii) = (n_cntrl-3) + l*(ii-1) + jj;
    end
    
end % end segment iteration loop

%% Assemble Column Vectors %%
    
iRowB_Eq = [ reshape(Dbnd0_iRow,[n_fix0*l 1]);
          reshape(Dbndf_iRow,[n_fixf*l 1]);
          reshape(Duhat_iRow,[3*n_seg 1]) ];
      
jColB_Eq = [ reshape(Dbnd0_jCol,[n_fix0*l 1]);
          reshape(Dbndf_jCol,[n_fixf*l 1]); 
          reshape(Duhat_jCol,[3*n_seg 1]) ];
    