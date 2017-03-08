function [zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt)
% function [zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt)
% 
% This function breaks up the column vector of design variables used to
% solve the collocation problem into several 3D matrices containing the
% state, control, and slack variables. If the collocation problem utilizes
% coast parameters then these are also output along with the number of
% design variables contained within a segment.
%
% INPUTS:
%    Z              design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    colt           structure containing collocation and optimization parameters
%
% OUTPUTS:
%    zis            design variable matrix (l x 1 x n_seg)
%    xis            state variable matrix (n_state x (N+1)/2 x n_seg)
%    uis            control variable matrix (n_cntrl x 1 x n_seg)
%    sis            slack variable matrix (n_slack x 1 x n_seg)
%    coast_times    vector of coast time parameters (n_coast x 1)
%    l              number of design variables per segment ((N+1)*n_state/2 + n_cntrl + n_slack)
%
% Originally Written by: R. Pritchett, 11/10/2016
% Last Update: R. Pritchett, 11/10/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N;
n_seg = colt.n_seg;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
n_coast = colt.n_coast;

% Remove coast and manifold prop times from design variable vector
coast_times = Z(1:n_coast); % will equal zero if n_coast = 0
Z = Z(n_coast+1:end); % states and controls at thrust arc nodes

% Define useful index numbers for state, control, and slack variables
l = (N+1)*n_state/2 + n_cntrl + n_slack ; % number of design variables per segment
u_start = 1; % initial index of control variables in zis
u_end = u_start+n_cntrl-1; % final index of control variables in zis
s_start = u_end+1; % initial index of slack variables in zis
s_end = s_start+n_slack-1; % final index of slack variables in zis

% Define logical indices 
index_state = true(1, l);
index_state(u_start:s_end) = false;
index_cntrl = false(1, l);
index_cntrl(u_start:u_end) = true;
index_slack = false(1, l);
index_slack(s_start:s_end) = true;

% Reshape Z into a 3D matrix 
zis = reshape(Z,[l,1,n_seg]);

% Break zis into state, control, and slack variables
xis = reshape(zis(index_state,:,:), [n_state (N+1)/2 n_seg]);
uis = repmat(zis(index_cntrl,:,:), [1 (N+1)/2 1]);
sis = repmat(zis(index_slack,:,:), [1 (N+1)/2 1]);
