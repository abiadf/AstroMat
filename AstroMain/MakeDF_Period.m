function [DF] = MakeDF_Period(Z,STMi,mult)
% function [DF] = MakeDF_Period(Z,STMi,mult)
%
% This function calculates the Jacobian matrix for multiple shooting 
% problems. For efficiency the Jacobian matrix is defined as a sparse
% matrix, thus this function calculates only the partials that are not 
% guaranteed to be nonzero values. However, some value in the output DF 
% vector may still be zero.
%
% INPUTS:
%    Z          design variable vector, (n*n_state x 1)
%    STMi       3D matrix of STM between relating each subsequent pair of
%               patch points
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    DF         nonzeros of the Jacobian matrix in column vector form
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
period = mult.period;
n_period = mult.n_period;
n_hypr0 = mult.n_hypr0;
mu = mult.mu;

% Convert column vector of design variables into matrix of state variables
x_ppt = reshape(Z(1:end-1),[n_state n]);

% Calculate patch point times using total time
T = Z(end);
t_ppt = linspace(0,T,n);
t_diff = repmat(diff(t_ppt),[n_state 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partials with Respect to Continuity Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Assemble STM Components %
%-------------------------------------------------------------------------%
 
% Reshape STMi into single column vector
DF_STM = reshape(STMi,[n_state*n_state*(n-1) 1]);

%-------------------------------------------------------------------------%
% Assemble Negative Identity Matrix Components %
%-------------------------------------------------------------------------%

% Preallocate storage matrix
neg_eye = zeros(n_state,n_state,n-1);

% Create 3D matrix of negative identity matrices
for ii = 1:n-1
    neg_eye_i = -eye(n_state);
    neg_eye(:,:,ii) = neg_eye_i;
end

% Reshape 3D matrix of identity matrices into single column vector
DF_neg_eye = reshape(neg_eye,[n_state*n_state*(n-1) 1]);

%-------------------------------------------------------------------------%
% Calculate and Assemble Partials with Respect to Time %
%-------------------------------------------------------------------------%
 
% Calculate partials with respect to time at every patch point
x_ppt_dt = EOM_CR3BP_vec(x_ppt,mu)./n; % necessary to divide by n because design variable is total time

% Reshape into single column vector
DF_dt = reshape(x_ppt_dt(:,2:end),[n_state*(n-1) 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partials with Respect to Boundary Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Calculate Partials with Respect to Periodicity Constraints %
%-------------------------------------------------------------------------%

% Partials of periodicity constraint with respect to patch point partials
DF_prdc0 = -ones(n_period,1);

% Partials of periodicity constraint with respect to propagated states at the final time
DF_prdcf_dx = reshape(STMi(period,:,end),[n_period*n_state 1]);

% Partials of periodicity constraint with respect to time variable
DF_prdcf_dt = x_ppt_dt(period,end);

% Combine partials of periodicity constraint into single column vector
DF_prdc = [DF_prdc0; DF_prdcf_dx; DF_prdcf_dt];

%-------------------------------------------------------------------------%
% Calculate Partials with Respect to Hyperplane Constraint %
%-------------------------------------------------------------------------%
DF_hypr0 = ones(n_hypr0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine all partials into single output vector
DF = [DF_STM; DF_neg_eye; DF_dt; DF_prdc; DF_hypr0];
