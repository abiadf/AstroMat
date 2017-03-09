function [t_bnd_new,maxe] = CEPadd(x_bnd,uis,t_bnd,colt)
% function [t_bnd_new,maxe] = CEPadd(x_bnd,t_bnd,n_seg,n_state,add_tol,mu)
% 
% This codes performs the segment addition portion of Control Explicit 
% Propagation (CEP) mesh refinement on a given set of boundary point states
% and their associated times.
%
% IMPORTANT VARIABLES:
%
% INPUTS:
%    x_bnd   Matrix of boundary point states (n_state x 1 x m)
%    t_bnd   Times at  boundary points (m x 1)
%    colt    structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    Tbnd_new   Times at new boundary points (m x 1)
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
add_tol = colt.add_tol;
options = colt.options;
ce = colt.ce;
mu = colt.mu;

% Calculate segment time intervals
t_seg = diff(t_bnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Segments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate Xend matrix
Xend = zeros(n_state,1,n_seg+1);
Xend(:,:,1) = x_bnd(:,:,1);

% Integrate to obtain states at each boundary points
ind = 1; % initialize iteration counter
for ii = 1:n_seg
    
    ind = ind+1; % iterate iteration counter
    
    % Set integration initial values
    x_bnd0 = x_bnd(:,:,ii);
    
    % Define constant control values for current segment
    cntrl = uis(:,1,ii);
    
    % Define integration timespan
    tspan = [0 t_seg(ii)];
    
    % Integrate
    [~,Xend_i] = ode113(@(t,in)EOM_CR3BP_LTCSI(t,in,cntrl,ce,mu),tspan,x_bnd0,options);
    Xend(:,:,ind) = Xend_i(end,:).';
    
end

%Obtain norm of the error at each boundary point
err = Xend - x_bnd;
err_norm = sqrt(sum(err.^2,1));

%Calculate max error for output to table
maxe = max(err_norm);

%Store indices of points where error is over tolerance
i_err = find(err_norm > add_tol);

%Calculate the number of new points that will be added
n_add = length(i_err);

%Preallocate new matrices for boundary point times
t_bnd_new = zeros(n_seg+1+n_add,1);

%Initialize vectors of indices for old and new points
i_o2n = 1:n_seg+1;
i_add = [i_err; 0]; % append by 0 to simplify for loop iteration

%Determine indices for old and new points
for jj = 1:n_add
    
    %Set index of current boundary point
    i_e = i_err(jj);
    
    %Iterate indices of old boundary points
    i_o2n(i_e:end) = i_o2n(i_e:end) + 1;
    
    %Iterate indices of new bounday points
    i_add(jj+1:end) = i_add(jj+1:end) + 1;
    
    %Interpolate times to obtain new boundary points
    t = t_bnd(i_e-1,1) + (t_bnd(i_e,1) - t_bnd(i_e-1,1)) ./ 2;
    t_bnd_new(i_add(jj),1) = t;

end

% %Remove appended value from i_add
% i_add = i_add(1:end-1);
    
%Add old boundary point states into Xbnd_new
for k = 1:n_seg+1
    t_bnd_new(i_o2n(k),1) = t_bnd(k,1); 
end
