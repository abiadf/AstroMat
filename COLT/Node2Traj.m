function [x_traj_store,t_traj_store,traj_TorC] = Node2Traj(x_bnd,uis,t_bnd,colt)
% function [x_traj,t_traj] = ppt2Traj(x_ppt,t_ppt,mult)
% 
% This function propagates between the boundary nodes of the solution to a
% collocation problem. It outputs 
%
% INPUTS:
%    x_bnd      nondimensional states at boundary nodes, (n_state x n)
%    t_ppt      nondimensional times at boundary nodes, (1 x n)
%    uis        3D matrix of control states for each segment (n_cntrl x (N+1)/2 x n)
%    colt       structure containing collocation parameters
%
% OUTPUTS:
%    x_traj_store     nondimensional states at integration time steps along entire trajectory
%    t_traj_store     nondimensional times at integration time steps along entire trajectory
%    traj_TorC        array of 'T' and 'C' characters that indicates 
%                     whether each segment is a thrust or a coast arc
%
% Written by R. Pritchett, 02/10/17
% Last Update: R. Pritchett, 02/10/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg = colt.n_seg; % initialize new segment value
options = colt.options;
ce = colt.ce;
mu = colt.mu;

% Calculate segment time intervals
t_seg = diff(t_bnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Segments %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate storage cell arrays
x_traj_store = cell(n_seg,1);
t_traj_store = cell(n_seg,1);
traj_TorC = cell(n_seg,1);

% Integrate to obtain states at each boundary points
for ii = 1:n_seg
    
    % Set integration initial values
    x_bnd0 = x_bnd(ii,:);
    
    % Define constant control values for current segment
    cntrl = uis(:,1,ii);
    
    % Calculate dimensional thrust magnitude for current segment
    Tmag_dim = uis(1,1,ii)*colt.m_ch*(1000*colt.l_ch)/(colt.t_ch)^2;
    
    % Define integration timespan
    tspan = [0 t_seg(ii)];
    
    % Integrate
    [t_traj,x_traj] = ode113(@(t,in)EOM_CR3BP_LTCSI(t,in,cntrl,ce,mu),tspan,x_bnd0,options);
    
    % Store numerical integration results
    x_traj_store{ii,1} = x_traj;
    t_traj_store{ii,1} = t_traj;

    % Record segment as thrust or coast arc
    tol = 1e-6; % if the thrust value is lower than this tolerance it is considered a coast arc
    if Tmag_dim <= tol
        traj_TorC{ii,1} = 'C';
    else
        traj_TorC{ii,1} = 'T';
    end 
    
end

