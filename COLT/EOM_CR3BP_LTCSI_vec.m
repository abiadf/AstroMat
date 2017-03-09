function [xis_dot] = EOM_CR3BP_LTCSI_vec(xis,uis,ce,mu)
% function [indot] = EOM_CR3BP_LTCSI_vec(xis,uis,mu)
%
% This function contains the equations of motion for a spacecraft equipped 
% with a constant specific impulse low-thrust engine in the circular 
% restricted three body problem. It is intended to be used when the input, 
% xis and uis are 2D or 3D matrices of state and control variables. 
%
% INPUTS:
%    xis    2D or 3D matrix of state values, each row must correspond to 
%           a different state, i.e. row 1 = x
%    uis    2D or 3D matrix of control values, each row must correspond to 
%           a different control variable, i.e. row 1 = Tmag
%    ce     effective exhaust velocity of CSI LT engine
%    mu     CR3BP mass ratio      
%
% OUTPUTS:
%    xis_dot  derivatives of input states
%
% Originally Written by: R. Pritchett, 09/19/2015
% Last Update: R. Pritchett, 09/29/2016
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

% Acceleration due to Thrust
aT1 = uis(1,:,:)./xis(7,:,:).*uis(2,:,:);
aT2 = uis(1,:,:)./xis(7,:,:).*uis(3,:,:);
aT3 = uis(1,:,:)./xis(7,:,:).*uis(4,:,:);

% Calculate EOM
xis_dot = [                 xis(4,:,:) ;
                            xis(5,:,:) ;
                            xis(6,:,:) ;
 2.*xis(5,:,:) + xis(1,:,:) + ad.*dht1 + ar.*rht1 + aT1;
-2.*xis(4,:,:) + xis(2,:,:) + ad.*dht2 + ar.*rht2 + aT2;
                             ad.*dht3 + ar.*rht3 + aT3;
                           -1.*uis(1,:,:)./ce ]; 
