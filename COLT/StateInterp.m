function [x] = StateInterp(t,C,t_pts,N)
% function [x] = StateInterp(t,C,t_pts,N)
% 
% This codes uses the converged solution of the existing mesh to
% interpolate the states x at a given time t. The existing mesh is given by
% a vector of boundary node times t_pts. The coefficients of the Nth degree
% polynomials interpolating between these boundary nodes are given by C 
%
% IMPORTANT VARIABLES:
%
% INPUTS:
%    t      dimensional time along trajectory at which to interpolate 
%    C      matrix of polynomial coefficients (l x (N+1) x n_seg)
%    ctol   Error tolerance for de Boor's method 
%    t_pts  vecor of boundary node dimensional times
%    N      degree of collocation polynomial
%
% OUTPUTS:
%    x      states at input time t
%
% See Grebow and Pavlak (2015) for details on the algorithm 
%cd

% Written by R. Pritchett, 7/17/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine which segment t is contained by and the times at its bounday
%nodes
lb_ind = find(t_pts<t);
lb_ind = lb_ind(end);
lb_t = t_pts(lb_ind);
ub_t = t_pts(lb_ind+1);
seg_ind = lb_ind;

%Calculate the corresponding nondimensional time
tau = 2/(ub_t-lb_t)*(t-lb_t)-1;

%Create B vector
B = (ones(1,N+1).*tau).^(0:N);

%Multiply B by C at the appropriate segment to obtain x
x = C(:,:,seg_ind)*B.';
