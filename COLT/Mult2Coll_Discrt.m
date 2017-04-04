function [Z0,t_var,t_bnd,xis_bnd,xis_plot,uis_plot] = Mult2Coll_Discrt(x_ppt,t_ppt,colt)
% function [Z0,t_var,t_bnd,xis_bnd,xis_plot,uis_plot] = Mult2Coll_Discrt(x_ppt,t_ppt,colt)
%
% This function discretizes an initial guess for a low-thrust transfer 
% trajectory into boundary and variable nodes such that the trajectory can 
% be passed to a collocation method. Specifically, this function
% discretizes revolutions about two different periodic orbits and combines
% these revolutions into a single initial guess.
%
% INPUTS:
%    colt     structure containing collocation and optimization parameters      
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

% Extract initial and final orbit IC
X_init = OrbIC(1:6);
T_init = OrbIC(7);
X_fin = OrbIC(8:13);
T_fin = OrbIC(14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Boundary and Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Boundary Nodes %%%

% Define boundary node times for each phase of the initial guess
t_bnd = t_ppt;

% Calculate segment timespans
t_seg = diff(t_bnd);

%%% Variable Nodes %%%

% Calculate nondimensional normalized time for each segment
[tau_nodes,~,~,~,~] = CollSetup(N,NodeSpace); % LG is interpolation method

% Preallocate t_var_vec matrix
t_var_vec = zeros((N+1)/2,1,n_seg);

for ii = 1:n_seg
    
    % Define time span for current segment
    t_bnd_i = [0 t_seg(ii)];
    n_seg_i = 1;

    % Compute all variable node dimensional times
    [t_var_i] = CollVarTimes(t_bnd_i,N,n_seg_i,tau_nodes,NodeSpace,0);
    
    % Store variable node times in 3D matrix
    t_var_vec(:,:,ii) = t_var_i; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate to Obtain Variable Node States %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate x_var_vec matrix
x_var_vec = zeros(6,(N+1)/2,n_seg);

% Calculate states at variable node times
for ii = 1:n_seg
    
    % Select variable node times for current segment
    t_var_i = t_var_vec(:,:,ii);
    
    % Select initial states for current segment
    x0_seg = x_ppt(ii,:);
    
    % Propagate to obtain boundary nodes states of initial orbit revolutions
    [x_var_i,~] = GslInteg_CR3BP(x0_seg,[0 t_var_i'],mu);
    
    % Store variable node states in 3D matrix
    x_var_vec(:,:,ii) = x_var_i(2:end,:)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble Discretized Coast Arc Outputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
uis_var = repmat(U,[1 1 n_seg]);
sis_var = repmat(S,[1 1 n_seg]);

% Add mass state to state values for extra periodic orbit revolutions
xis_bnd = [x_ppt(:,1:6) m0.*ones(n_seg+1,1)]';
xis_var = cat(1,x_var_vec, m0.*ones(1,(N+1)/2,n_seg));

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
    xis_var_i = xis_var(:,:,ii);
    uis_var_i = uis_var(:,:,ii);
    sis_var_i = sis_var(:,:,ii);
    
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
%% Shift Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_bnd_vec = repmat(reshape(t_bnd,[1 1 n_seg+1]),[(N+1)/2 1 1]);
t_var_shift = t_var_vec + t_bnd_vec(:,:,1:n_seg);
t_var = reshape(t_var_shift,[n_seg*(N+1)/2 1]);
