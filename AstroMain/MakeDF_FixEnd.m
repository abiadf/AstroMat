function [DF] = MakeDF_FixEnd(Z,STMi,mult)
% function [DF] = MakeDF_FixEnd(Z,STMi,mult)
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
mu = mult.mu;
i_hypr0 = mult.i_hypr0; % defines which component of the initial state is constrained to a hyperplane 
x_hypr0 = mult.x_hypr0; % defines the value(s) of the initial state component(s) that define hyperplane 
n_hypr0 = mult.n_hypr0; % defines the number of initial states constrained to a hyperplane
i_hyprf = mult.i_hyprf; % defines which component of the final state is constrained to a hyperplane 
x_hyprf = mult.x_hyprf; % defines the value(s) of the final state component(s) that define hyperplane 
n_hyprf = mult.n_hyprf; % defines the number of final states constrained to a hyperplane

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
% Calculate Partials with Respect to Hyperplane Constraints %
%-------------------------------------------------------------------------%

DF_hypr0 = ones(n_hypr0,1);
DF_hyprf = ones(n_hyprf,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine all partials into single output vector
DF = [DF_STM; DF_neg_eye; DF_dt; DF_hypr0; DF_hyprf];
