function [DF] = MakeDF_Num(Z,Zl,Fl,mult)
% function [DF] = MakeDF_Num(Z,Zl,Fl,mult)
% 
% This function calculates the Jacobian matrix for multiple shooting
% problems using forward or central difference methods. This function is
% intended to be used to check the results of the numerically integrated
% Jacobian only and not to be the primary method of computing the Jacobian 
% in the multiple shooting algorithm.
%
% INPUTS:
%    Z          design variable vector, (n_state*n+1 x 1)
%    Zl         number of design variables
%    Fl         number of constraints
%    mult       structure containing multiple shooting parameters
%
% OUTPUTS:
%    DF         Jacobian matrix
%
% Written by R. Pritchett, 02/09/17
% Last Update: R. Pritchett, 02/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from coll stucture
DiffType = mult.DiffType;

% Pre-allocate matrices
DF = zeros(Fl,Zl); % DF matrix
pert = zeros(Zl,1); % perturbation matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute DF Using Forward Step or Central Difference Differentiation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch DiffType
    
    %---------------------------------------------------------------------%
    % Compute DF Using Forward Step Differentiation %
    %---------------------------------------------------------------------%
    case 'Forward'
        
    % Set step size for forward step differentiation
    stpsz = sqrt(eps);    

        for ii = 1:Zl % loop through design variables

%             display(ii)
            
            % Simultaneously perturb 1 state of 1 node
            pert(ii,1)=stpsz; 

            % Add/Subtract perturbation to variable node states
            Z_pls = Z+pert;
    
            % Compute constraints with plus perturbation        
            [F_pls,~] = MakeF(Z_pls,mult);

            % Compute constraints without perturbation
            [F,~] = MakeF(Z,mult);

            % Apply forward difference differentiation
            DF(:,ii)=(F_pls-F)./(stpsz);
    
            pert=0*pert; % reset perturbation matrix
            
        end
        
    %---------------------------------------------------------------------%
    % Compute DF Using Central Difference Differentiation %
    %---------------------------------------------------------------------%    
    case 'Central'
        
        % Set step size for central difference approximation
        stpsz = eps^(1/3); 
    
        for ii = 1:Zl % loop through design variables
            
            % Simultaneously perturb 1 state of 1 node
            pert(ii,1)=stpsz; 

            % Add/Subtract perturbation to variable node states
            Z_pls = Z+pert;
            Z_min = Z-pert;
    
            % Compute constraints with plus perturbation        
            [F_pls,~] = MakeF(Z_pls,mult);

            % Compute constraints without perturbation
            [F_min,~] = MakeF(Z_min,mult);

            % Apply central difference differentiation
            DF(:,ii)=(F_pls-F_min)./(2*stpsz);
    
            pert=0*pert; % reset perturbation matrix
            
        end
end