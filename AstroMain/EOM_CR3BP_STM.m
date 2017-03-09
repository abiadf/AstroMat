function [indot] = EOM_CR3BP_STM(t,in,mu)
% function [indot] = EOM_CR3BP_STM(t,in,mu)
% 
% This function contains the ballistic equations of motion for a particle 
% in the CR3BP expressed in the rotating frame as well as the first order 
% differential equations for the state transition matrix. It is intended 
% to be used with one of Matlab's built-in ode integrators. 
%
% INPUTS:
%    t      integration time step, required for use with ode integrators 
%    in     input states [42x1]
%    mu     CR3BP mass ratio      
%
% OUTPUTS:
%    indot  derivatives of input states [42x1]
%
% Originally Written by: R. Pritchett, 09/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Breakup input vector
x = in(1);
y = in(2);
z = in(3);
xdot = in(4);
ydot = in(5);
zdot = in(6);
phi = reshape(in(7:42),[6 6]).';

%Calculate d and r
d = sqrt((x+mu)^2+y^2+z^2);
r = sqrt((x-1+mu)^2+y^2+z^2);

%Integrate Equations of Motion
rdot(1,1) = xdot;
rdot(1,2) = ydot;
rdot(1,3) = zdot;
rdot(1,4) = 2*ydot+x-(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3;
rdot(1,5) = -2*xdot+y-(1-mu)*y/d^3-mu*y/r^3;
rdot(1,6) = -(1-mu)*z/(d^3)-(mu*z)/(r^3);

%Integrate State Transition Matrix
Uxx = 1-(1-mu)/(d^3)-mu/(r^3)+(3*(1-mu)*(x+mu)^2)/(d^5)+(3*mu*(x-1+mu)^2)/(r^5);
Uxy = (3*(1-mu)*(x+mu)*y)/(d^5)+(3*mu*(x-1+mu)*y)/(r^5);
Uxz = (3*(1-mu)*(x+mu)*z)/(d^5)+(3*mu*(x-1+mu)*z)/(r^5);
Uyx = Uxy;
Uyy = 1-(1-mu)/(d^3)-mu/(r^3)+(3*(1-mu)*y^2)/(d^5)+(3*mu*y^2)/(r^5);
Uyz = (3*(1-mu)*y*z)/(d^5)+(3*mu*y*z)/(r^5);
Uzx = Uxz;
Uzy = Uyz;
Uzz = -(1-mu)/(d^3)-mu/(r^3)+(3*(1-mu)*(z^2))/(d^5)+(3*mu*z^2)/(r^5);
row1 = [0 0 0 1 0 0]; 
row2 = [0 0 0 0 1 0];
row3 = [0 0 0 0 0 1];
row4 = [Uxx Uxy Uxz 0 2 0];
row5 = [Uyx Uyy Uyz -2 0 0];
row6 = [Uzx Uzy Uzz 0 0 0];

%Assemble A matrix
A = [row1;row2;row3;row4;row5;row6];

%Compute phidot
phidot = A*phi;

%Break phidot down into rows for output
phidot = reshape(phidot.',[1 36]);
indot = [rdot phidot]'; 



