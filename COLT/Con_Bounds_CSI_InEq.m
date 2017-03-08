function [FB_InEq] = Con_Bounds_CSI_InEq(uis,colt)
% function [FB_InEq] = Con_Bounds_CSI_InEq(uis,colt)
% 
% This function computes the boundary inequality constraints for the 
% constraint vector in the low-thrust transfer problem with a CSI engine 
% This constraint vector includes constraints on the minimum and maximum
% values of the thrust magnitude.
%
% INPUTS:
%    uis            control variable matrix (n_cntrl x (N+1)/2 x n_seg)
%    colt           structure containing collocation and optimization parameters
%
% OUTPUTS:
%    FB_InEq        boundary inequality constraint vector
%
% Originally Written by: R. Pritchett, 02/01/2016
% Last Update: R. Pritchett, 02/02/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg = colt.n_seg;
Tmax = colt.Tmax; 

% ConB Inequality 1: Thrust constrained via slack variables to be in range: 0<T<Tmax
FB_InEq(1,:,:) = -uis(1,1,:);
FB_InEq(2,:,:) = uis(1,1,:) - Tmax;
FB_InEq = reshape(FB_InEq,[2*n_seg,1]);
