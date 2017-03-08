function [ZPlot] = Z2Plot(Z,x_bnd,t_bnd,t_var,colt)
% function [ZPlot] = Z2Plot(Z,x_bnd,t_bnd,colt)
% 
% This function extracts data from the design variable vector Z and
% processes it such that it can be easily plotted. The state, control, 
% and slack variable quantities from the design variable vector are broken
% into separate variables. Additionally, if any coast parameters are
% included in the design variables these values are used to propagate the
% appropriate coast arcs.
%
% INPUTS:
%    Z        design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    x_bnd    matrix of state values at boundary nodes (n_state x n_seg+1)
%    t_bnd    non-normalized times at boundary nodes (n_seg x 1)
%    t_var    non-normalized times at variable nodes (n_seg*(N+1)/2)
%    colt     structure containing collocation and optimization parameters
%
% OUTPUTS:
%    ZPlot    structure containing vectors and matrices for plotting
%
% Written by R. Pritchett, 10/04/16
% Last Update: R. Pritchett, 10/04/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
OrbIC = colt.OrbIC;
phi0 = colt.phi0;
options = colt.options;
N = colt.N;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
n_coast = colt.n_coast;
PorM = colt.PorM;
ce = colt.ce;
ce_dim = colt.ce_dim;
Isp = colt.Isp;
Isp_dim = colt.Isp_dim;
mu = colt.mu;

% Calculate new number of segments and save to structure
n_seg = length(t_bnd)-1;

% If there are coast parameters extract times from vector
if n_coast == 2 % Only periodic orbit coast times
    tau1 = Z(1); % coast time along initial periodic orbit 
    tau2 = Z(2); % coast time along final periodic orbit  
elseif n_coast == 4 % Periodic orbit and manifold coast times
    tau1 = Z(1); % coast time along initial periodic orbit 
    tau2 = Z(2); % coast time along final periodic orbit
    alph1 = Z(3); % coast time for manifold 1 
    alph2 = Z(4); % coast time for manifold 2 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Breakup Design Variable Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert column vector of design variables to 3D matrix
[zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);

% Reshape output for plotting (No 3D matrices)
x_bnd_plot = reshape(x_bnd,[7 n_seg+1])'; % states at boundary nodes
x_var_plot = reshape(xis,[7 n_seg*(N+1)/2])'; % states at variable nodes
u_var = reshape(uis(2:4,:,:),[3 n_seg*(N+1)/2])'; % thrust unit vector
T_var = reshape(uis(1,:,:),[1 n_seg*(N+1)/2])'; % thrust magnitude
uT_var = repmat(T_var,[1 3]).*u_var; % full thrust vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate Between Boundary Nodes to Obtain Thrust and Coast Segments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through nodes and propagate
[x_traj,t_traj,traj_TorC] = Node2Traj(x_bnd_plot,uis,t_bnd,colt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Dimensional State and Control Variable Values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate dimensional values of state variables
pos_bnd_dim_plot = x_bnd_plot(:,1:3).*colt.l_ch; % [km]
vel_bnd_dim_plot = x_bnd_plot(:,4:6).*(colt.l_ch/colt.t_ch); % [km/s]
mass_bnd_dim_plot = x_bnd_plot(:,7).*colt.m_ch; % [kg]
pos_var_dim_plot = x_var_plot(:,1:3).*colt.l_ch; % [km]
vel_var_dim_plot = x_var_plot(:,4:6).*(colt.l_ch/colt.t_ch); % [km/s]
mass_var_dim_plot = x_var_plot(:,7).*colt.m_ch; % [kg] 

% Calculate dimensional values of control variables
T_var_dim = T_var.*colt.m_ch*(1000*colt.l_ch)/(colt.t_ch)^2; % [Newtons]
uT_var_dim = uT_var.*colt.m_ch*(1000*colt.l_ch)/(colt.t_ch)^2; % [Newtons]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reshape Control Variable Vectors to Enable Histogram Plotting Style %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modify vector of variable node times
t_var_plot = reshape(repmat(t_var',[2,1]),[1 2*length(t_var)]);
t_var_hist = [t_bnd(1) t_var_plot(2:end-1) t_bnd(end)];

% Modify vector of thrust magnitude
T_var_plot = reshape(repmat(T_var_dim(1:end-1)',[2,1]),[1 2*(length(T_var_dim)-1)]);
T_var_hist = [T_var_plot(1) T_var_plot T_var_plot(end)];

% Modify thrust pointing vectors one component at a time
uT_var_hist = zeros(3,2*length(uT_var_dim)); 
for ii = 1:3
    uT_var_plot = reshape(repmat(uT_var_dim(1:end-1,ii)',[2,1]),[1 2*(length(uT_var_dim)-1)]);
    uT_var_hist(ii,:) = [uT_var_plot(:,1) uT_var_plot uT_var_plot(:,end)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagate Any Relevant Coast Arcs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If used plot coast arcs along periodic orbits and/or manifold arcs 
if n_coast == 2
    
    % Propagate to obtain new periodic orbit coast arcs
    [~,X_tau1] = ode113(@(t,in)EOM_STM(t,in,mu),[0 tau1],[X_init phi0],options);
    [~,X_tau2] = ode113(@(t,in)EOM_STM(t,in,mu),[0 tau2],[X_fin phi0],options);

elseif n_coast == 4
    
    % Propagate to obtain new periodic orbit coast arcs
    [~,X_tau1] = ode113(@(t,in)EOM_CR3BP(t,in,mu),[0 tau1],OrbIC(1:6),options);
    [~,X_tau2] = ode113(@(t,in)EOM_CR3BP(t,in,mu),[0 tau2],OrbIC(8:13),options);

    % Propagate manifolds arcs to obtain initial and final thrust points
    [X_alph1] = ManiProp(tau1,alph1,'initial','unstable',PorM{1},colt);
    [X_alph2] = ManiProp(tau2,alph2,'final','stable',PorM{2},colt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Results in Structure for Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define structure for output
ZPlot = struct;

% Save results in structure for output to plotting script
ZPlot.x_traj = x_traj;
ZPlot.t_traj = t_traj;
ZPlot.traj_TorC = traj_TorC;
ZPlot.x_bnd_plot = x_bnd_plot;
ZPlot.x_var_plot = x_var_plot;
ZPlot.u_var = u_var;
ZPlot.T_var = T_var;
ZPlot.uT_var = uT_var;
ZPlot.pos_bnd_dim_plot = pos_bnd_dim_plot;
ZPlot.vel_bnd_dim_plot = vel_bnd_dim_plot;
ZPlot.mass_bnd_dim_plot = mass_bnd_dim_plot;
ZPlot.pos_var_dim_plot = pos_var_dim_plot;
ZPlot.vel_var_dim_plot = vel_var_dim_plot;
ZPlot.mass_var_dim_plot = mass_var_dim_plot;
ZPlot.T_var_dim = T_var_dim;
ZPlot.uT_var_dim = uT_var_dim;
ZPlot.t_var_hist = t_var_hist;
ZPlot.T_var_hist = T_var_hist;
ZPlot.uT_var_hist = uT_var_hist;

% If coast parameters are used save the corresponding trajectory
if n_coast == 2
    ZPlot.X_tau1 = X_tau1;
    ZPlot.X_tau2 = X_tau2;
    
elseif n_coast == 4
    ZPlot.X_tau1 = X_tau1;
    ZPlot.X_tau2 = X_tau2;
    ZPlot.X_alph1 = X_alph1;
    ZPlot.X_alph2 = X_alph2;
end
    

