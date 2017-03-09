function [indot] = EOM_CR3BP_LTCSI(t,in,cntrl,ce,mu)
% function [indot] = EOM_CR3BP_LTCSI(t,in,cntrl,mu)
% 
% This function contains the equations of motion for a spacecraft equipped 
% with a constant specific impulse low-thrust engine in the circular 
% restricted three body problem. The four control variables are thrust 
% magnitude and the three components of thrust direction. These are assumed 
% to be constant over the duration of the propagation. This function is 
% intended to be used with Matlab's built-in ode integrators. 
%
% INPUTS:
%    t        integration time step, must be included to enable use with Matlab ode integrators 
%    in       input states
%    cntrl    vector of control variables; values are constant throughout integration
%    mu       CR3BP mass ratio      
%
% OUTPUTS:
%    indot    derivatives of input states
%
% Originally Written by: R. Pritchett, 02/01/2017
% Last Update: R. Pritchett, 02/01/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breakup state vector
x = in(1);
y = in(2);
z = in(3);
xdot = in(4);
ydot = in(5);
zdot = in(6);
m = in(7);

% Breakup vector of constant control parameters
T = cntrl(1); % thrust magnitude
u = cntrl(2:4); % thrust direction (3x1) vector

%% Position, Velocity, and Mass Process Equations %%

% Calculate d and r
d = sqrt((x+mu)^2+y^2+z^2);
r = sqrt((x-1+mu)^2+y^2+z^2);

% Integrate Equations of Motion
indot(1:3,1) = in(4:6);
indot(4,1) = 2*ydot+x-(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3 + T/m*u(1);
indot(5,1) = -2*xdot+y-(1-mu)*y/d^3-mu*y/r^3 + T/m*u(2);
indot(6,1) = -(1-mu)*z/(d^3)-(mu*z)/(r^3) + T/m*u(3);
indot(7,1) = -T./ce;
