function [FC,STMi,x_fin] = Con_Cont(x_ppt,t_ppt,mult)
% function [FC] = Con_Cont(x_ppt,t_ppt,mult)
% 
% This function calculates the continuity constraints for multiple shooting
% problems.
%
% INPUTS:
%    x_ppt      nondimensional states at patch points, (n_state x n)
%    t_ppt      nondimensional times at patch points, (1 x n)
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    FC         continuity constraint vector
%    STMi       3D matrix of STM between relating each subsequent pair of
%               patch points
%    x_fin      states at final patch point obtained from propagation,
%               needed for computation of boundary constraints
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
phi0 = mult.phi0;
mu = mult.mu;

% Preallocate storage matrices
FC = zeros(n-1,n_state);
STMi = zeros(n_state,n_state,n-1);

% Calculate timespan between patch points
t_ppt_diff = diff(t_ppt);

% Propagate between patch points and compute constraint violations
for ii = 1:n-1
    
    % Define segment initial conditions
    x0_seg = [x_ppt(:,ii)', phi0];
    
    % Define segment timespan
    tspan = [0 t_ppt_diff(ii)];
    
    % Propagate
    [x_seg,~] = GslInteg_STM(x0_seg,tspan,mu);
    
    % Define states and STM at final time
    xf_seg = x_seg(end,1:n_state); % states at final time
    phif_seg = x_seg(end,n_state+1:end); % STM at final time
    
    % Calculate constraint
    FC(ii,:) = xf_seg - x_ppt(:,ii+1)';
    
    % Save STM for assembling Jacobian 
    phif_seg = reshape(phif_seg,[n_state n_state])';
    STMi(:,:,ii) = phif_seg;
    
end

% Reshape FC into single column vector
FC = reshape(FC',[n_state*(n-1) 1]);

% Output states at final patch point obtained from propagation
x_fin = xf_seg;
    
    


