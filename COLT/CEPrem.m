function [t_bnd_new,maxe] = CEPrem(x_bnd,uis,t_bnd,colt)
% function [t_bnd_new,maxe] = CEPrem(x_bnd,Tbnd,n_seg,n_state,rem_tol,mu)
% 
% This codes performs the segment removal portion of Control Explicit Propagation (CEP) mesh refinement 
% on a given set of boundary point states and their associated times.
%
% IMPORTANT VARIABLES:
%
% INPUTS:
%    x_bnd   Matrix of boundary point states (n_state x 1 x m)
%    t_bnd   Times at  boundary point (m x 1)
%    colt      structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    t_bnd_new  Times at new boundary points (m x 1)
%    maxe      maximum error at boundary points 
%
% See Grebow and Pavlak (2015) for details on the algorithm
%
% Written by R. Pritchett, 10/20/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg = colt.n_seg; % initialize new segment value
n_state = colt.n_state;
rem_tol = colt.rem_tol;
options = colt.options;
ce = colt.ce;
mu = colt.mu;

% Calculate segment time intervals
t_seg = diff(t_bnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Segments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if even or odd number of segments
rem_chk = rem(n_seg,2);

% Adjust for loop end condition according to even or odd result 
if rem_chk == 0 % if even
    for_end = n_seg-1;
else % if odd
    for_end = n_seg-2; 
end

% Calculate indices of odd points
i_odd = 1:2:n_seg+1;

% Calculate number of odd points
n_odd = length(i_odd);

% Preallocate Xend_odd matrix
Xend_odd = zeros(n_state,1,n_odd);
Xend_odd(:,:,1) = x_bnd(:,:,1);

% Propagate to obtain states at odd boundary points
ind = 1;
for jj = 1:2:for_end
    
    ind = ind+1;
    
    % Set integration initial values
    x_bnd0 = x_bnd(:,:,jj);
    
    % Define constant control values for current segment
    cntrl = uis(:,1,jj);
    
    % Define integration timespan
    tspan = [0 sum(t_seg(jj:jj+1,1))];
    
    % Integrate
    [~,Xend] = ode113(@(t,in)EOM_CR3BP_LTCSI(t,in,cntrl,ce,mu),tspan,x_bnd0,options);
    Xend_odd(:,:,ind) = Xend(end,:).';
    
end
    
% Create matrix of states at odd boundary points
x_bnd_odd = x_bnd(:,:,i_odd);

% Obtain norm of the error at each odd boundary point
err = Xend_odd - x_bnd_odd;
err_norm = sqrt(sum(err.^2,1));

% Calculate max error for output to table
maxe = max(err_norm);

% Store indices of points where error is under tolerance
i_err = 2*(find(err_norm < rem_tol))-1;

% Calculate the number of points that will be removed
n_rem = length(i_err);

% Create array of indices of points to remove
i_rem = i_err-1;

% Create array of indices of points to be kept
i_keep = 1:n_seg+1;
for k = 1:n_rem
    ir = i_rem(k);
    i_keep = i_keep(i_keep~=ir);
end

% Create matrix of times for points to be kept
t_bnd_new = t_bnd(i_keep,1);

