function [indot] = EOM_CR3BP(t,in,mu)
% function [indot] = EOM_STM(t,in,mu)
% 
% This function contains the equations of motion of the circular restricted
% three body problem. It is intended to be used with Matlab's 
% built-in ode integrators 
%
% INPUTS:
%    t      integration time step, must be included to enable use with Matlab ode integrators 
%    in     input states
%    mu     CR3BP mass ratio      
%
% OUTPUTS:
%    indot  derivatives of input states
%
% Originally Written by: R. Pritchett, 01/15/2016
% Last Update: R. Pritchett, 01/15/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Breakup input vector
x = in(1);
y = in(2);
z = in(3);
xdot = in(4);
ydot = in(5);
zdot = in(6);

%Calculate d and r
d = sqrt((x+mu)^2+y^2+z^2);
r = sqrt((x-1+mu)^2+y^2+z^2);

%Integrate Equations of Motion
indot(1,1) = xdot;
indot(2,1) = ydot;
indot(3,1) = zdot;
indot(4,1) = 2*ydot+x-(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3;
indot(5,1) = -2*xdot+y-(1-mu)*y/d^3-mu*y/r^3;
indot(6,1) = -(1-mu)*z/(d^3)-(mu*z)/(r^3);



