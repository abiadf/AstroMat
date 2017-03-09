function [Z0,x_ppt,t_ppt,x_plot] = L4Short_Discretize(T,mult)
% function [Z0,x_ppt,t_ppt,x_plot] = L4Short_Discretize(T,mult)
%
% This function discretizes a short period L4 orbit by integration timesteps, 
% boundary nodes, and variable nodes so that the trajectory can 
% be passed as an initial guess to a collocation method. 
%
% INPUTS:
%    T        orbital period
%    mult     structure containing multiple shooting parameters parameters      
%
% OUTPUTS:
%    Z0       column vector containing initial guess for a multiple 
%             shooting method.
%    x_ppt    matrix of patch point states 
%    t_ppt    vector of patch point times
%    x_plot   states at patch points collected in matrix form for simple plotting 
%
% Originally Written by: R. Pritchett, 02/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
options = mult.options;
L = mult.L;
mu = mult.mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Boundary and Variable Node Times %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Boundary Nodes %%%

% Define boundary node times
t_ppt = linspace(0,T,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Parameters of Linear Variational Equations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate 2nd Partials of the Pseudo-Potential Function at L4
Uxx = 3/4;
Uyy = 9/4;
Uxy = -(3*sqrt(3)/2)*(mu-0.5); % choose negative because calculating at L4 (positive if calculating at L5)

% Calculate s1 and s2 frequencies
g = 1 - 27*mu*(1-mu);
s1 = imag(sqrt(0.5*(-1+sqrt(g))));
s2 = imag(sqrt(0.5*(-1-sqrt(g))));

% Define beta coefficients based on desired motion
b1 = 0; % set to 1 to obtain long period motion
b2 = 0; % set to 1 to obtain long period motion
b3 = 0.05; % set to 1 to obtain short period motion
b4 = 0.05; % set to 1 to obtain short period motion

% Calculate alpha coefficients
a1 = (-2*b2*s1 - Uxy*b1)/(Uxx+s2^2);
a2 = (-2*b1*s1 + Uxy*b2)/(Uxx+s2^2);
a3 = (-2*b4*s2 - Uxy*b3)/(Uxx+s2^2);
a4 = (-2*b3*s2 + Uxy*b4)/(Uxx+s2^2);

% Define phasing angle (radians)
phi = 3*pi/4; % with this phasing the trajectory begins at y = y_L4
% phi = 7*pi/4; % with this phasing the trajectory begins at y = y_L5
% phi = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Variable Node States and Assemble Design Variable Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solutions of the Variational Equations near the Equilateral Points
xi_var = a1.*cos(s1.*t_ppt + phi) + a2.*sin(s1.*t_ppt + phi) + a3.*cos(s2.*t_ppt + phi) + a4.*sin(s2.*t_ppt + phi);
eta_var = b1.*cos(s1.*t_ppt + phi) + b2.*sin(s1.*t_ppt + phi) + b3.*cos(s2.*t_ppt + phi) + b4.*sin(s2.*t_ppt + phi);
xi_dot_var = -a1*s1.*cos(s1.*t_ppt + phi) + a2*s1.*sin(s1.*t_ppt + phi) - a3*s2.*cos(s2.*t_ppt + phi) + a4*s2.*sin(s2.*t_ppt + phi);
eta_dot_var = -b1*s1.*cos(s1.*t_ppt + phi) + b2*s1.*sin(s1.*t_ppt + phi) - b3*s2.*cos(s2.*t_ppt + phi) + b4*s2.*sin(s2.*t_ppt + phi);

% Add location of L4 Lagrange point to variations
x_var = xi_var' + L(4,1).*ones(n,1);
y_var = eta_var' + L(4,2).*ones(n,1);

% Define zero vectors for z and zdot states since motion is planar
z_var = zeros(n,1);
z_dot_var = zeros(n,1);

% Assemble matrix of states at variable nodes
x_ppt = [x_var y_var z_var xi_dot_var' eta_dot_var' z_dot_var];

% Reshape states into column vector of design variables
Z0 = reshape(x_ppt',[n_state*n 1]);

% Add total trajectory time to design variable vector
Z0 = [Z0;T];

% Reshape (if necessary) x_ppt into form amenable to plotting
x_plot = x_ppt;