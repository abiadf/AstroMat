function [FB] = Con_Bounds_CSI(uis,sis,x0,xf,colt)
% function [FB] = Con_Bounds_CSI(uis,sis,x0,xf,coast_times,colt)
% 
% This function computes the boundary constraints for the constraint vector
% in the low-thrust transfer problem with a VSI engine and coast time 
% parameters This constraint vector includes constraints on continuity, 
% intial mass, thrust pointing vector magnitude, the sign of the thrust 
% magnitude, and maximum power.
%
% INPUTS:
%    uis            control variable matrix (n_cntrl x (N+1)/2 x n_seg)
%    sis            slack variable matrix (n_slack x (N+1)/2 x n_seg)
%    x0             initial states of each segment (n_state x 1 x n_seg)
%    xf             final states of each segment (n_state x 1 x n_seg)
%    colt           structure containing collocation and optimization parameters
%
% OUTPUTS:
%    FB     boundary constraint vector
%
% Originally Written by: R. Pritchett, 06/07/2016
% Last Update: R. Pritchett, 09/29/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N;
n_seg = colt.n_seg;
n_state = colt.n_state;
n_cntrl = colt.n_cntrl;
n_slack = colt.n_slack;
Tmax = colt.Tmax; 
x0_des = colt.x0_des;
x0_des_ind = colt.x0_des_ind;
xf_des = colt.xf_des;
xf_des_ind = colt.xf_des_ind;

% Ensure x0_des and xf_des are column vectors
x0_des = x0_des(:);
xf_des = xf_des(:);

% ConB 1: Endpoints constrained to be fixed for specified components
Fbnd0 = x0(1:n_state,1,1) - x0_des; % calculate difference between initial endpoint and desired initial endpoint
Fbnd0 = Fbnd0(x0_des_ind); % keep only initial endpoint components that are specified as fixed
Fbndf = xf(1:n_state,1,end) - xf_des; % calculate difference between final endpoint and desired initial endpoint
Fbndf = Fbndf(xf_des_ind); % keep only final endpoint components that are specified as fixed

% ConB 2: Pointing vector constrained to be unit magnitude
uhat = reshape(uis(n_cntrl-2:n_cntrl,1,:),[3 n_seg]);
Fumag = (sum(uhat.^2)).^(1/2) - ones(1,n_seg);
Fumag = Fumag.';

% ConB 3: Thrust constrained via slack variables to be in range: 0<T<Tmax
Fslack(1,:,:) = uis(1,1,:) - sis(1,1,:).^2; % bound T to be greater than 0
Fslack(2,:,:) = Tmax - uis(1,1,:) - sis(2,1,:).^2; % bound T to be less than Tmax
Fslack = reshape(Fslack,[n_slack*n_seg,1]);

% Assemble boundary constraints into single vector
FB = [Fbnd0; Fbndf; Fumag; Fslack];