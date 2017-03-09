function [xis_dot] = EOM_CR3BP_vec(xis,mu)
% function [xis_dot] = EOM_CR3BP_vec(xis,mu)
% 
% This function contains the equations of motion of the circular restricted
% three body problem and is intended to be used when the input, xis, is a 
% 2D or 3D matrix of state variables. 
%
% INPUTS:
%    xis      2D or 3D matrix of state values, each row must correspond to 
%             a different state, i.e. row 1 = x 
%    mu       CR3BP mass ratio      
%
% OUTPUTS:
%    xis_dot  derivatives of input states
%
% Originally Written by: R. Pritchett, 02/06/2017
% Last Update: R. Pritchett, 02/06/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the magnitude and components of r
rmg=sqrt((xis(1,:,:)+mu-1).^2+xis(2,:,:).^2+xis(3,:,:).^2);
rht1=(xis(1,:,:)+mu-1)./rmg;
rht2=xis(2,:,:)./rmg;
rht3=xis(3,:,:)./rmg;

% Calculate the magnitude and components of d
dmg=sqrt((xis(1,:,:)+mu).^2+xis(2,:,:).^2+xis(3,:,:).^2);
dht1=(xis(1,:,:)+mu)./dmg;
dht2=xis(2,:,:)./dmg;
dht3=xis(3,:,:)./dmg;

% Acceleration due to P1
ad = -(1-mu)./dmg.^2;

% Acceleration due to P2
ar = -mu./rmg.^2;

% Calculate EOM
xis_dot = [                 xis(4,:,:) ;
                            xis(5,:,:) ;
                            xis(6,:,:) ;
 2*xis(5,:,:) + xis(1,:,:) + ad.*dht1 + ar.*rht1;
-2*xis(4,:,:) + xis(2,:,:) + ad.*dht2 + ar.*rht2;
                             ad.*dht3 + ar.*rht3 ]; 
