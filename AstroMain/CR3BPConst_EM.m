% CR3BPConst_EM.m
%
% This script calculates various constants used in the Circular 
% Restricted Three Body Problem (CR3BP) as applied to the Earth-Moon 
% system. It is intended to be run at the beginning of a main script so 
% that all necessary variables will be available in the workspace. 
% Note that the parameter values used in this script were obtained from 
% JPL's Monte software.
%
% Originally Written by: R. Pritchett, 07/24/2015
% Last Update: R. Pritchett, 02/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define Constants %% 

% Define gravitational constant and acceleration due to gravity on the
% surface of the Earth
G = 6.67384e-20; % gravitational constant [km^3/(kg*s^2)]
g0 = 9.8196; % acceleration due to gravity [m/s^2]

% Define radi and semi-major axes of primary bodies
r_E = 6371; % radius of the Earth, [km]
r_M = 1737.53; % radius of the Moon, [km]
a_M = 382918.1893584794; % semi-major axis of Moon orbit from Monte at epoch 1-jan-2015 et, [km]

% Define GM values of primary bodies. Notes, these are the values used in JPL's Monte Software
mu_earth = 3.986004362333397e5; % Earth gravitational parameter [km^3/s^2]
mu_moon = 4.902800076227743e3; % Moon gravitational parameter [km^3/s^2]
mu = mu_moon/(mu_earth+mu_moon); % mu defined for CR3BP [nondimensional]
M = (mu_earth+mu_moon)/G;

%% Characteristic Quantities %%

l_ch = a_M; %Characteristic length
m_ch = M; %Characteristic mass
t_ch = sqrt((l_ch^3)/(mu_earth+mu_moon)); %Characteristic time
v_ch = l_ch/t_ch; %Characteristic velocity

%% Calculate Libration Points %%

% Run function that calculates libration point positions
L = libration_pos_EM(mu);

%% Define Initial State Transition Matrices %%

% Create initial phi matrix, used to compute the STM for ballistic motion
phi0 = reshape(eye(6,6),[1 36]);

% Create initial psi matrix, used to compute the STM for low-thrust motion
psi0 = reshape(eye(14).',[1 196]);
