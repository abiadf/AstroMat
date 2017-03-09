function [Z0,t_var,t_bnd,xis_plot,uis_plot] = L2LNearHeteroDiscrt_InEq(X,colt)
% function [Z0,t_var,t_bnd,xis_plot,uis_plot] = L2LNearHeteroDiscrt_InEq(X,colt)
%
% This function discretizes a low-thrust transfer trajectory into boundary 
% nodes, and variable nodes so that the trajectory can be passed as an 
% initial guess to a collocation method. It can be applied to low-thrust 
% trajectories obtained via indirect optimization, or near heteroclinic
% transfers. Note that this version of the discretization function does 
% NOT include slack variables because it assumes the initial guess will be
% passed to a function, such as fmincon, capable of handling inequality
% constraints. This is why the function name is appended with _InEq.
%
% INPUTS:
%    X        vector of time parameter values,(4x1)
%    colt     structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    Z0        column vector containing initial guess for collocation
%              method
%    t_var     row vector of variable node times, 1 x n_seg*(N+1)/2 
%    t_bnd     row vector of boundary node times, 1 x n_seg+1 
%    xis_plot  matrix of variable node states composed in a form that is
%              conducive to plotting
%    uis_plot  matrix of variable node control states composed in a form 
%              that is conducive to plotting
%
% Originally Written by: R. Pritchett, 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
N = colt.N;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
phi0 = colt.phi0;
m0 = colt.m0;
options = colt.options;
mu = colt.mu;
OrbIC = colt.OrbIC;
l_ch = colt.l_ch;

% Extract initial and final orbit IC
X_init = OrbIC(1:6);
T_init = OrbIC(7);
X_fin = OrbIC(8:13);
T_fin = OrbIC(14);

% Breakup design variable vector from indirect optimization solution
tau0 = X(1); % coast time on orbit 1
alph0 = X(2); % propagation time for manifold 1 
alphf = X(3); % propagation time for manifold 2 
tauf = X(4); % coast time on orbit 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Boundary and Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define number of segments for each phase of the initial guess
n_seg_revi = 5; % number of segments for initial orbit rev
n_seg_P1 = 10; % number of segments for Phase 1, initial orbit coast
n_seg_P2 = 50; % number of segments for Phase 2, unstable manifold
n_seg_P3 = 50; % number of segments for Phase 3, stable manifold
n_seg_P4 = 10; % number of segments for Phase 4, final orbit coast
n_seg_revf = 5; % number of segments for final orbit rev
n_seg = n_seg_revi+n_seg_P1+n_seg_P2+n_seg_P3+n_seg_P4+n_seg_revf; % total number of segments, initially

% Define propagation times for each phase of transfer initial guess
t_revi = T_init;
t_P1 = tau0; % transfer time for phase 1
t_P2 = alph0; % transfer time for phase 2
t_P3 = abs(alphf); % transfer time for phase 3
t_P4 = T_fin - tauf; % transfer time for phase 4
t_revf = T_fin;

%%% Boundary Nodes %%%

% Define boundary node times
t_bnd_revi = linspace(0,t_revi,n_seg_revi+1);
t_bnd_P1 = linspace(0,t_P1,n_seg_P1+1);
t_bnd_P2 = linspace(0,t_P2,n_seg_P2+1);
t_bnd_P3 = linspace(0,-t_P3,n_seg_P3+1); % Note: will integrate in reverse time
t_bnd_P4 = linspace(0,-t_P4,n_seg_P4+1); % Note: will integrate in reverse time
t_bnd_revf = linspace(0,t_revf,n_seg_revf+1);

% Define times in total transfer time for each phase
TF_revi = t_revi; % total time after initial orbit rev
TF_P1 = TF_revi + t_P1; % total time after initial orbit coast
TF_P2 = TF_P1 + t_P2; % total time after unstable manifold coast phase
TF_P3 = TF_P2 + t_P3; % total time after stable manifold coast phase
TF_P4 = TF_P3 + t_P4; % total time after final orbit coast phase
TF_revf = TF_P4 + t_revf; % total time after final orbit rev

%%% Variable Nodes %%%

% Calculate nondimensional normalized time for each segment
[tau_nodes,~,~,~,~] = CollSetup(N,NodeSpace); % LG is interpolation method

% Compute all variable node dimensional times
[t_var_revi] = CollVarTimes(t_bnd_revi,N,n_seg_revi,tau_nodes,NodeSpace,0);
[t_var_P1] = CollVarTimes(t_bnd_P1,N,n_seg_P1,tau_nodes,NodeSpace,0);
[t_var_P2] = CollVarTimes(t_bnd_P2,N,n_seg_P2,tau_nodes,NodeSpace,0);
[t_var_P3] = CollVarTimes(t_bnd_P3,N,n_seg_P3,tau_nodes,NodeSpace,0);
[t_var_P4] = CollVarTimes(t_bnd_P4,N,n_seg_P4,tau_nodes,NodeSpace,0);
[t_var_revf] = CollVarTimes(t_bnd_revf,N,n_seg_revf,tau_nodes,NodeSpace,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate Initial and Final Orbit Extra Revolutions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propagate to obtain variable nodes of initial orbit extra revolution
[~,X_revi_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_revi],[X_init phi0],options);
X_revi_var = X_revi_var(2:end,:); % t=0 not variable node so remove states

% Propagate to obtain variable nodes of final orbit extra revolution
[~,X_revf_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_revf],[X_fin phi0],options);
X_revf_var = X_revf_var(2:end,:); % t=0 not variable node so remove states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate Periodic Orbit Coast Arcs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propagate to obtain variable nodes of unstable orbit coast arc
[~,X_tau1_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_P1],[X_init phi0],options);
X_tau1_var = X_tau1_var(2:end,:); % t=0 not variable node so remove states

% Propagate to obtain variable nodes of stable orbit coast arc; Recall this is propagated in REVERSE time
[~,X_tau2_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_P4],[X_fin phi0],options);
X_tau2_var = X_tau2_var(2:end,:); % t=0 not variable node so remove states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Manifold Step-Off Points %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define monodromy matrix for initial and final orbits
mono_orb1 = reshape(X_revi_var(end,7:42),[6 6])';
mono_orb2 = reshape(X_revf_var(end,7:42),[6 6])';

% Calculate eigenvalues and eigenvectors for monodromy matrices
[V1_mono,~] = eig(mono_orb1); % note the transpose of mono_orb1 is the true monodromy matrix
[V2_mono,~] = eig(mono_orb2); % note the transpose of mono_orb2 is the true monodromy matrix

% Propagate along orbits 1 and 2 for coast times tau1 and tau2
[~,X_tau1] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0 tau0],[X_init phi0],options);
[~,X_tau2] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0 tauf],[X_fin phi0],options);

% Set step-off distance along the eigenspace
d_dim = 50; % [km]
d = d_dim/l_ch;

% Create STM's at step-off points
Phi_tau1 = reshape(X_tau1(end,7:42),[6 6])';
Phi_tau2 = reshape(X_tau2(end,7:42),[6 6])';    
    
% Obtain all 6 eigenvectors at step-off points
evec1 = Phi_tau1*V1_mono; 
evec2 = Phi_tau2*V2_mono;

% Define individual eigenvectors corresponding to stable and unstable eigenvalues
vU = evec1(:,1); % only want unstable eigenspace from Orbit 1
vS = evec2(:,2); % only want stable eigenspace from Orbit 2

% Calculate stable and unstable eigenspaces
Wu = vU/norm(vU);
Ws = vS/norm(vS);
    
% Calculate states of step-off point on half-manifold
xUplus = X_tau1(end,1:6)' + d.*Wu;
% xUminus = X_tau1(end,1:6)' - d.*Wu;
xSplus = X_tau2(end,1:6)' + d.*Ws;
% xSminus = X_tau2(end,1:6)' - d.*Ws;

% Define which half-manifold step-off occurs on
xm1 = xUplus';
% xm1 = xUminus';
xm2 = xSplus';
% xm2 = xSminus';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate Manifolds %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propagate to obtain variable nodes of unstable manifold
[~,X_alph1_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_P2],[xm1 phi0],options);
X_alph1_var = X_alph1_var(2:end,:); % t=0 not variable node so remove states

% Propagate to obtain variable nodes of unstable manifold; Recall this is propagated in REVERSE time
[~,X_alph2_var] = ode113(@(t,in)EOM_CR3BP_STM(t,in,mu),[0; t_var_P3],[xm2 phi0],options);
X_alph2_var = X_alph2_var(2:end,:); % t=0 not variable node so remove states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Discretized Coast Arc Outputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define final mass equal to initial mass, no mass loss in initial guess
mf = m0;

% Define control values along coast arcs
T = 0; % zero thrust corresponds to coast phase
uhat = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3)]; % thrust pointing vector is always unit magnitude
U = [T uhat]';

% Set control and slack values for extra periodic orbit revolutions
uis_var_revi = repmat(U,[1 1 n_seg_revi]);
uis_var_revf = repmat(U,[1 1 n_seg_revf]);

% Set periodic orbit coast arc control and slack values
uis_var_P1 = repmat(U,[1 1 n_seg_P1]);
uis_var_P4 = repmat(U,[1 1 n_seg_P4]);

% Set manifold coast arc control and slack values
uis_var_P2 = repmat(U,[1 1 n_seg_P2]);
uis_var_P3 = repmat(U,[1 1 n_seg_P3]);

% Add mass state to state values for extra periodic orbit revolutions
X_Urev_var = [X_revi_var(:,1:6) m0.*ones((N+1)/2*n_seg_revi,1)];
X_Srev_var = [X_revf_var(:,1:6) mf.*ones((N+1)/2*n_seg_revf,1)];

% Add mass state to periodic orbit coast arc state values
X_Ucst_var = [X_tau1_var(:,1:6) m0.*ones((N+1)/2*n_seg_P1,1)];
X_Scst_var = [X_tau2_var(:,1:6) mf.*ones((N+1)/2*n_seg_P4,1)];

% Add mass state to manifold coast arc state values
X_Uman_var = [X_alph1_var(:,1:6) m0.*ones((N+1)/2*n_seg_P2,1)];
X_Sman_var = [X_alph2_var(:,1:6) mf.*ones((N+1)/2*n_seg_P3,1)];

% Reshape vector of state values at variable nodes along initial periodic orbit extra revolution
xis_var_revi = zeros(n_state,(N+1)/2,n_seg_revi);
ind = 1;
for ii = 1:n_seg_revi
    xis_var_revi(:,:,ii) = X_Urev_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along final periodic orbit extra revolution
xis_var_revf = zeros(n_state,(N+1)/2,n_seg_revf);
ind = 1;
for ii = 1:n_seg_revf
    xis_var_revf(:,:,ii) = X_Srev_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along initial periodic orbit coast arc
xis_var_P1 = zeros(n_state,(N+1)/2,n_seg_P1);
ind = 1;
for ii = 1:n_seg_P1
    xis_var_P1(:,:,ii) = X_Ucst_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along final periodic orbit coast arc
xis_var_P4 = zeros(n_state,(N+1)/2,n_seg_P4);
X_Scst_var_flip = flip(X_Scst_var,1);
ind = 1;
for ii = 1:n_seg_P4
    xis_var_P4(:,:,ii) = X_Scst_var_flip(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along unstable manifold
xis_var_P2 = zeros(n_state,(N+1)/2,n_seg_P2);
ind = 1;
for ii = 1:n_seg_P2
    xis_var_P2(:,:,ii) = X_Uman_var(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

% Reshape vector of state values at variable nodes along stable manifold
xis_var_P3 = zeros(n_state,(N+1)/2,n_seg_P3);
X_Sman_var_flip = flip(X_Sman_var,1);
ind = 1;
for ii = 1:n_seg_P3
    xis_var_P3(:,:,ii) = X_Sman_var_flip(ind:ind+(N-1)/2,:)';
    ind = ind + (N+1)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Full Discretized Halo to Halo Transfer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine state matrices into single matrix
xis_var_TF = cat(3,xis_var_revi,xis_var_P1,xis_var_P2,xis_var_P3,xis_var_P4,xis_var_revf);

% Combine control matrices into single matrix
uis_var_TF = cat(3,uis_var_revi,uis_var_P1,uis_var_P2,uis_var_P3,uis_var_P4,uis_var_revf);

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
    
    % Formulate single column for Z_mat
    x_col = reshape(xis_var_i,[n_state*(N+1)/2 1]);
    Z_mat(:,ii) = [uis_var_i; x_col];
    
    % Assemble variables into plotting friendly form
    xis_plot(:,ind:ind+(N-1)/2) = xis_var_i;
    uis_plot(:,ii) = uis_var_i;
    ind  = ind + (N+1)/2;
    
end

% Reshape Z_mat into column vector
Z0 = reshape(Z_mat,[l*n_seg 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shift Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total t_bnd_vector
t_bnd = [t_bnd_revi TF_revi+t_bnd_P1(2:end) TF_P1+t_bnd_P2(2:end) TF_P2+(-t_bnd_P3(2:end)) TF_P3+(-t_bnd_P4(2:end)) TF_P4+t_bnd_revf(2:end)]';

% Total t_var vector
t_var = [t_var_revi; TF_revi+t_var_P1; TF_P1+t_var_P2; TF_P2+(-t_var_P3); TF_P3+(-t_var_P4); TF_P4+t_var_revf];
