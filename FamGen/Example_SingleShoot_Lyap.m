% Example_Lyapunov_EM.m
%
% This script uses collocation and direct optimization to compute a
% Lyapunov orbit 
%
% Originally Written by: R. Pritchett, 02/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc;

%% Set Constants %%

%Call script that calculates useful CR3BP constants
CR3BPConst_EM

% Define collocation and optimization parameter structure
sshoot = struct;

% Save CR3BP constant parameters in structure 
sshoot.mu = mu; % gravitational 
sshoot.l_ch = l_ch; % characteristic length
sshoot.t_ch = t_ch; % characteristic time
sshoot.L = L; % libration point position coordinates 
sshoot.phi0 = phi0; % initial phi matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Numerical Integration Options %%%
RelativeTol = 1e-12; % relative integration tolerance
AbsoluteTol = 1e-12; % absolute integration tolerance
sshoot.options = odeset('RelTol',RelativeTol,'AbsTol',AbsoluteTol);
sshoot.options_event = odeset('Events',@xcross,'RelTol',1e-12,'AbsTol',1e-12);

%%% Shooting Inputs %%%
sshoot.newt_tol = 1e-12; % tolerance for Newton's Method
sshoot.atten_tol = 1e-4; % tolerance for attenuation factor
sshoot.atten = 1; % attenuation factor for Newton's method step size 1/atten
sshoot.DiffType = 'Forward'; % finite difference method for numerical jacobian
sshoot.n_state = 6; % number of state variables
sshoot.n_slack = 0; % number of slack variables

%%% Targeter Case Inputs %%%
sshoot.ContCase = 'NoCont'; % target single orbit, no continuation
% sshoot.ContCase = 'PArc'; % target orbit for use in pseudo-arclength arclength continuation scheme
% sshoot.ContCase = 'NatParam'; % target orbit for use in natural parameter continuation
sshoot.OrbCase = 'Lyapunov'; % target Lyapunov or other planar orbits symmetric about the x-axis
% sshoot.OrbCase = 'Halo'; % target halo orbits
% sshoot.OrbCase = 'Vertical'; % target vertical orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set of initial conditions from Tom Pavlak
load EM_L2_Lyap_250_orbits

% Define initial conditions
orb_num = 50;
X0 = u_0_save(orb_num,:);
T0 = t_f_save(orb_num);
T0_half = T0/2;
tspan = [0 T0];

% Integrate initial guess
[x_init,t_init] = GslInteg_CR3BP(X0,tspan,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal

plot(X0(1),X0(2),'*k');
plot(x_init(:,1),x_init(:,2),':r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve Single Shooting Problem %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run multiple shooting algorithm
[X,T_half,DF_fin] = PeriodicTarget(X0,T0_half,sshoot);

% Save or load single shooting results
% save('SSLyap','X','T_half')
% load SSLyap.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Final Result Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integrate converged result guess
[x_fin,t_fin] = GslInteg_CR3BP(X,[0 2*T_half],mu);

% Plot converged result
plot(x_fin(:,1),x_fin(:,2),'b');
hold off
