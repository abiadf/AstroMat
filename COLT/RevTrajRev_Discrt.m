function [Z0,t_var,t_bnd,xis_bnd,xis_plot,uis_plot] = RevTrajRev_Discrt(traj_guess,colt)
% function [Z0,t_var,t_bnd,xis_bnd,xis_plot,uis_plot] = RevTrajRev_Discrt(traj_guess,colt)
%
% This function discretizes an initial guess for a low-thrust transfer 
% trajectory into boundary and variable nodes such that the trajectory can 
% be passed to a collocation method. Specifically, this function
% discretizes revolutions about two different periodic orbits and combines
% these revolutions with an intermediate trajectory arc into a single
% initial guess.
%
% INPUTS:
%    colt        structure containing collocation and optimization parameters
%    traj_guess  structure containing initial guess information 
%
% OUTPUTS:
%    Z0          column vector containing initial guess for collocation
%                method
%    t_var       column vector of variable node times
%    t_bnd       column vector of boundary node times
%    xis_bnd     matrix of boundary node states
%    xis_plot    matrix of variable node states configured for simple plotting 
%    uis_plot    matrix of variable node controls configured for simple
%                plotting
%
% Originally Written by: R. Pritchett, 09/23/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from traj_guess stucture
TrajType = traj_guess.TrajType;
Z_ConstOrInterp = traj_guess.Z_ConstOrInterp;

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
N = colt.N;
n_seg = colt.n_seg;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
m0 = colt.m0;
mu = colt.mu;
OrbIC = colt.OrbIC;
n_rev = colt.n_rev;
n_seg_revi = colt.n_seg_revi;
n_seg_revf = colt.n_seg_revf;
n_seg_trans = colt.n_seg_trans;
% n_seg = n_seg_revi + n_seg_trans + n_seg_revf;

% Extract initial and final orbit IC
X_init = OrbIC(1:6);
T_init = OrbIC(7);
X_fin = OrbIC(8:13);
T_fin = OrbIC(14);

% Extract number of revs to propagate for initial and final orbits
n_revi = n_rev(1);
n_revf = n_rev(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretize Revolutions about Initial and Final Orbit %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Calculate Boundary and Variable Node Times
%-------------------------------------------------------------------------%

% Define propagation times for each phase of the initial guess
t_revi = n_revi*T_init;
t_revf = n_revf*T_fin;

%%% Boundary Nodes %%%

% Define boundary node times for each phase of the initial guess
t_bnd_revi = linspace(0,t_revi,n_seg_revi+1);
t_bnd_revf = linspace(0,t_revf,n_seg_revf+1);

%%% Variable Nodes %%%

% Calculate nondimensional normalized time for each segment
[tau_nodes,~,~,~,~] = CollSetup(N,NodeSpace); % LG is interpolation method

% Compute all variable node dimensional times
[t_var_revi] = CollVarTimes(t_bnd_revi,N,n_seg_revi,tau_nodes,NodeSpace,0);
[t_var_revf] = CollVarTimes(t_bnd_revf,N,n_seg_revf,tau_nodes,NodeSpace,0);

%-------------------------------------------------------------------------%
% Propagate to Obtain Boundary and Variable Node States %
%-------------------------------------------------------------------------%

% Propagate to obtain boundary nodes states of initial orbit revolutions
[X_revi_bnd,~] = GslInteg_CR3BP(X_init,t_bnd_revi,mu);

% Propagate to obtain boundary node states of final orbit revolutions
[X_revf_bnd,~] = GslInteg_CR3BP(X_fin,t_bnd_revf,mu);

% Propagate to obtain variable node states of initial orbit revolutions
[X_revi_var,~] = GslInteg_CR3BP(X_init,[0 t_var_revi'],mu);
X_revi_var = X_revi_var(2:end,:); % t=0 not variable node so remove states

% Propagate to obtain variable node states of final orbit revolutions
[X_revf_var,~] = GslInteg_CR3BP(X_fin,[0 t_var_revf'],mu);
X_revf_var = X_revf_var(2:end,:); % t=0 not variable node so remove states

%-------------------------------------------------------------------------%
% Assemble Matrices of State, Control, and Slack Variables %
%-------------------------------------------------------------------------%

% Define final mass equal to initial mass, no mass loss in initial guess
mf = m0;

% Define control values along coast arcs
T = 0; % zero thrust corresponds to coast phase
uhat = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]; % thrust pointing vector is always unit magnitude
U = [T uhat]';

% Define slack variables along coast arcs
beta_Tlow = sqrt(0.1); % chosen such that slack variable constraint initially equals zero
beta_Thi = sqrt(0.1); % chosen such that slack variable constraint initially equals zero
S = [beta_Tlow beta_Thi]';

% Set control and slack values for extra periodic orbit revolutions
uis_var_revi = repmat(U,[1 1 n_seg_revi]);
uis_var_revf = repmat(U,[1 1 n_seg_revf]);
sis_var_revi = repmat(S,[1 1 n_seg_revi]);
sis_var_revf = repmat(S,[1 1 n_seg_revf]);

% Add mass state to state values for extra periodic orbit revolutions
X_revi_bnd = [X_revi_bnd(:,1:6) m0.*ones(n_seg_revi+1,1)];
X_revf_bnd = [X_revf_bnd(:,1:6) mf.*ones(n_seg_revf+1,1)];
X_revi_var = [X_revi_var(:,1:6) m0.*ones((N+1)/2*n_seg_revi,1)];
X_revf_var = [X_revf_var(:,1:6) mf.*ones((N+1)/2*n_seg_revf,1)];

% Reshape vectors of state values at boundary nodes
xis_bnd_revi = reshape(X_revi_bnd',[n_state n_seg_revi+1]);
xis_bnd_revf = reshape(X_revf_bnd',[n_state n_seg_revf+1]);

% Reshape vector of state values at variable nodes along initial periodic orbit extra revolution
xis_var_revi = zeros(n_state,(N+1)/2,n_seg_revi);
ind = 1;
for ii = 1:n_seg_revi
    xis_var_revi(:,:,ii) = X_revi_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along final periodic orbit extra revolution
xis_var_revf = zeros(n_state,(N+1)/2,n_seg_revf);
ind = 1;
for ii = 1:n_seg_revf
    xis_var_revf(:,:,ii) = X_revf_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretize Middle Trajectory Arc %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch TrajType
    
    case 'FromColl' % input trajectory arc is in a collocation format
        
        switch Z_ConstOrInterp
            
            case 'Const_Z' % Design variable vector remains the same
                
                % Extract necessary parameters from traj_guess stucture
                t_bnd_trans = traj_guess.t_bnd_old';
                t_var_trans = traj_guess.t_var_old;
                Z_trans = traj_guess.Z_trans;
                
                %-------------------------------------------------------------------------%
                % Convert Z into 3D Matrices and Rename Variables
                %-------------------------------------------------------------------------%

                % Convert column vector of design variables to 3D matrices
                colt.n_seg = colt.n_seg_trans; % change n_seg in colt to use Z23D function for Z_trans
                [~,xis_var_trans,uis_var_trans,sis_var_trans,~,~] = Z23D(Z_trans,colt);
                uis_var_trans = uis_var_trans(:,1,:); % only need first column since all columns are identical
                sis_var_trans = sis_var_trans(:,1,:); % only need first column since all columns are identical
                
                % Compute constant collocation matrices necessary to compute constraints
                [collmat] = OptSetup(t_bnd_trans,colt);

                % Run constraint function to calculate states at boundary nodes
                [~,x0,xf,~] = MakeF_LT(Z_trans,collmat.t_seg,collmat.t_seg_d,collmat,colt);
                colt.n_seg = n_seg; % revert n_seg in colt back to correct value

                %Create vector of boundary nodes 
                xis_bnd_trans = reshape(cat(3,x0,xf(:,:,end)),[n_state n_seg_trans+1]);
                
            case 'Interp_New_Z' % New design variable vector is interpolated 
        
                % Extract necessary parameters from traj_guess stucture
                n_seg_old = traj_guess.n_seg_old;
                t_bnd_old = traj_guess.t_bnd_old;
                t_var_old = traj_guess.t_var_old;
                Z_trans = traj_guess.Z_trans;
                C_trans = traj_guess.C_trans;
                tf_trans = traj_guess.tf_trans;

                %-------------------------------------------------------------------------%
                % Calculate Boundary Node Times and Setup for Interpolation
                %-------------------------------------------------------------------------%

                % Define boundary node times for trajectory arc
                t_bnd_trans = linspace(0,tf_trans,n_seg_trans+1);

                % Convert column vector of design variables to 3D matrices
                colt.n_seg = n_seg_old; % change n_seg in colt to use Z23D function for Z_trans
                [~,~,uis,sis,~,~] = Z23D(Z_trans,colt);
                colt.n_seg = n_seg; % revert n_seg in colt back to correct value

                %-------------------------------------------------------------------------%
                % Interpolate to Obtain State at Variable Nodes
                %-------------------------------------------------------------------------%

                % Interpolate to obtain new variable node states
                [Z_trans,t_var_trans] = ZInterp(t_bnd_old,t_bnd_trans,t_var_old,tau_nodes,...
                    C_trans,uis,sis,n_seg_old,n_seg_trans,colt); % Note: colt.n_seg can remain unchanged here because it is not used in ZInterp

                % Convert column vector of design variables to 3D matrices
                colt.n_seg = colt.n_seg_trans; % change n_seg in colt to use Z23D function for Z_trans
                [~,xis_var_trans,uis_var_trans,sis_var_trans,~,~] = Z23D(Z_trans,colt);
                uis_var_trans = uis_var_trans(:,1,:); % only need first column since all columns are identical
                sis_var_trans = sis_var_trans(:,1,:); % only need first column since all columns are identical

                %-------------------------------------------------------------------------%
                % Run Constraint Function to Calculate States at Boundary Nodes
                %-------------------------------------------------------------------------%

                % Compute constant collocation matrices necessary to compute constraints
                [collmat] = OptSetup(t_bnd_trans,colt);

                % Run constraint function to calculate states at boundary nodes
                [~,x0,xf,~] = MakeF_LT(Z_trans,collmat.t_seg,collmat.t_seg_d,collmat,colt);
                colt.n_seg = n_seg; % revert n_seg in colt back to correct value

                %Create vector of boundary nodes 
                xis_bnd_trans = reshape(cat(3,x0,xf(:,:,end)),[n_state n_seg_trans+1]);
                
        end % end Z_ConstOrInterp switch statement 
        
    case 'FromMult' %input trajectory arc is in a multiple shooting format
        
        % Extract necessary parameters from traj_guess stucture
        x_ppt = traj_guess.x_ppt;
        t_ppt = traj_guess.t_ppt;

        %-------------------------------------------------------------------------%
        % Calculate Boundary and Variable Node Times and States
        %-------------------------------------------------------------------------%

        %%% Boundary Nodes %%%
        
        % Define boundary node times for each phase of the initial guess
        t_bnd_trans = t_ppt';

        % Calculate segment timespans
        t_seg_trans = diff(t_bnd_trans);

        %%% Variable Nodes %%%
 
        % Preallocate t_var_vec matrix
        t_var_trans_vec = zeros((N+1)/2,1,n_seg_trans);
        x_var_vec = zeros(6,(N+1)/2,n_seg_trans);

        % Loop through segment to calculate variable node states and times
        for ii = 1:n_seg_trans

            % Define time span for current segment
            t_bnd_i = [0 t_seg_trans(ii)];
            n_seg_trans_i = 1;

            % Compute all variable node dimensional times
            [t_var_i] = CollVarTimes(t_bnd_i,N,n_seg_trans_i,tau_nodes,NodeSpace,0);

            % Store variable node times in 3D matrix
            t_var_trans_vec(:,:,ii) = t_var_i; 
            
            % Select initial states for current segment
            x0_seg = x_ppt(ii,:);

            % Propagate to obtain boundary nodes states of initial orbit revolutions
            [x_var_i,~] = GslInteg_CR3BP(x0_seg,[0 t_var_i'],mu);

            % Store variable node states in 3D matrix
            x_var_vec(:,:,ii) = x_var_i(2:end,:)';

        end
        
        %-------------------------------------------------------------------------%
        % Assemble Matrices of State, Control, and Slack Variables %
        %-------------------------------------------------------------------------%
        
        % Set control and slack values for extra periodic orbit revolutions
        uis_var_trans = repmat(U,[1 1 n_seg_trans]);
        sis_var_trans = repmat(S,[1 1 n_seg_trans]);

        % Add mass state to state values for extra periodic orbit revolutions
        xis_bnd_trans = [x_ppt(:,1:6) m0.*ones(n_seg_trans+1,1)]';
        xis_var_trans = cat(1,x_var_vec, m0.*ones(1,(N+1)/2,n_seg_trans));
        
        %-------------------------------------------------------------------------%
        % Reshape and Shift t_var_trans Vector %
        %-------------------------------------------------------------------------%
        
        t_bnd_trans_vec = repmat(reshape(t_bnd_trans,[1 1 n_seg_trans+1]),[(N+1)/2 1 1]);
        t_var_trans_shift = t_var_trans_vec + t_bnd_trans_vec(:,:,1:n_seg_trans);
        t_var_trans = reshape(t_var_trans_shift,[n_seg_trans*(N+1)/2 1]);

end % end TrajType switch statement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Full Discretized Transfer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine state matrices for boundary nodes into single matrix
xis_bnd = cat(2,xis_bnd_revi,xis_bnd_trans(:,2:end),xis_bnd_revf(:,2:end));

% Combine state matrices into single matrix
xis_var_TF = cat(3,xis_var_revi,xis_var_trans,xis_var_revf);

% Combine control matrices into single matrix
uis_var_TF = cat(3,uis_var_revi,uis_var_trans,uis_var_revf);

% Combine slack variable matrices into single matrix
sis_var_TF = cat(3,sis_var_revi,sis_var_trans,sis_var_revf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Design Variable Vector and Plotting Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = (N+1)*n_state/2 + n_cntrl + n_slack; % total number of design variable per segment
Z_mat = zeros(l, n_seg ); % preallocate matrix of design variable values
xis_plot = zeros(n_state, n_seg*(N+1)/2 ); % preallocate plotting friendly matrix of state variable values
uis_plot = zeros(n_cntrl, n_seg ); % preallocate plotting friendly matrix of state variable values
ind = 1; % initialize index counter for xis_plot and uis_plot

for ii = 1:n_seg
    
    % Define state, control, and slack variable for current segment
    xis_var_i = xis_var_TF(:,:,ii);
    uis_var_i = uis_var_TF(:,:,ii);
    sis_var_i = sis_var_TF(:,:,ii);
    
    % Formulate single column for Z_mat
    x_col = reshape(xis_var_i,[n_state*(N+1)/2 1]);
    Z_mat(:,ii) = [uis_var_i; sis_var_i; x_col];
    
    % Assemble variables into plotting friendly form
    xis_plot(:,ind:ind+(N-1)/2) = xis_var_i;
    uis_plot(:,ii) = uis_var_i;
    ind  = ind + (N+1)/2;
    
end

% Reshape Z_mat into column vector
Z0 = reshape(Z_mat,[l*n_seg 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shift Boundary and Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define time in terms of total transfer time for each phase of the initial guess
TF_revi = t_revi; % total time after initial orbit revs
TF_trans = TF_revi + t_bnd_trans(end);
TF_revf = TF_trans + t_revf; % total time after final orbit revs

% Total t_bnd_vector
t_bnd = [t_bnd_revi TF_revi+t_bnd_trans(2:end) TF_trans+t_bnd_revf(2:end)]';

% Total t_var vector
t_var = [t_var_revi; TF_revi+t_var_trans; TF_trans+t_var_revf];
