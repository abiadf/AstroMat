function [C] = Jacobi_Calc(x0,mu)
% function [C] = Jacobi_Calc(x0,mu)
% 
% This function calculates the Jacobi constant values corresponding to a
% matrix of n six states.
%
% INPUTS:
%    x0     matrix of n six states, (n x 6)
%    mu     CR3BP mass ratio      
%
% OUTPUTS:
%    C      vector of Jacobi constant values, (n x 1)
%
% Originally Written by: R. Pritchett, 02/23/2017
% Last Update: R. Pritchett, 02/23/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate nondimensional magnitude from primary and secondar bodies
dmg = sqrt((x0(:,1)+mu).^2 + x0(:,2).^2 + x0(:,3).^2);
rmg = sqrt((x0(:,1)-1+mu).^2 + x0(:,2).^2 + x0(:,3).^2);

% Calculate pseudopotential
Ustr = 0.5.*(x0(:,1).^2 + x0(:,2).^2) + (1-mu)./dmg + mu./rmg;

% Calculate velocity magnitude
vmg = x0(:,4).^2 + x0(:,5).^2 + x0(:,6).^2;

% Calculate Jacobi constant
C = 2.*Ustr - vmg;