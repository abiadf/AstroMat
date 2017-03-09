function [x_traj,t_traj] = ppt2Traj(x_ppt,t_ppt,mult)
% function [x_traj,t_traj] = ppt2Traj(x_ppt,t_ppt,mult)
% 
% This function propagates between patch points on a trajectory and 
% concatenates the resulting states into a single matrix that contains 
% states along the entire trajectory. The array of times resulting from 
% each propagation are treated in the same way  
%
% INPUTS:
%    x_ppt      nondimensional states at patch points, (n_state x n)
%    t_ppt      nondimensional times at patch points, (1 x n)
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    x_traj     nondimensional states at integration time steps along entire trajectory
%    t_traj     nondimensional times at integration time steps along entire trajectory
%
% Written by R. Pritchett, 02/10/17
% Last Update: R. Pritchett, 02/10/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
mu = mult.mu;

% Calculate timespan between patch points
t_ppt_diff = diff(t_ppt);

% Initialize x_traj and t_traj vectors
x_traj = [];
t_traj = [];

% Propagate between patch points and compute constraint violations
for ii = 1:n-1
    
    % Define segment initial conditions
    x0_seg = x_ppt(:,ii)';
    
    % Define segment timespan
    tspan = [0 t_ppt_diff(ii)];
    
    % Propagate
    [x_seg,t_seg] = GslInteg_CR3BP(x0_seg,tspan,mu);
    
    % Add results to total trajectory matrix
    x_traj = [x_traj; x_seg];
    t_traj = [t_traj; t_seg];
    
end