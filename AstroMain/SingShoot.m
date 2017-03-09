function [X,T_half,DF_fin,sshoot] = SingShoot(X0,T0_half,sshoot)
% function [X,T_half,DF_fin,sshoot] = SingShoot(X0,T0_half,sshoot)
% 
% Given an initial design variable vector, half period time, and a 
% structure of parameters this function solves a single shooting problem 
% for natural motion in the CR3BP. Specifically, it is meant to target 
% symmetric orbits in the CR3BP and allows different orbit cases and 
% continuation methods to be selected, these include: Lyapunov, Halo, 
% and Vertical Orbits.
%
% INPUTS:
%    Z0         initial design variable vector (n_state+1 x 1)
%    sshoot     structure containing single shooting parameters
%
% OUTPUTS:
%    Z          final design variable vector (n_state*n+1 x 1)
%    DF_fin     DF matrix used for final update needed for use in Pseudo-Arclength Continuation method
%    sshoot     structure containing single shooting parameters
%
% Written by R. Pritchett, 02/13/17
% Last Update: R. Pritchett, 02/13/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from sshoot stucture
OrbCase = sshoot.OrbCase;
ContCase = sshoot.ContCase; 

switch ContCase
    
    case 'NoCont'
    
        switch OrbCase

            case 'Lyapunov'

                % Run single shooting algorithm for Lyapunov and other planar
                % orbits symmetric about the x-axis
                [X,T_half,DF_fin,sshoot] = SingShoot_Lyap(X0,T0_half,sshoot);

            case 'Halo'
                
                % Run single shooting algorithm for Halo orbits and other
                % 3D orbits symmetric about the xz-plane
                [X,T_half,DF_fin,sshoot] = SingShoot_Halo(X0,T0_half,sshoot);

            case 'Vertical'
                               
        end
        
    case 'PArc'

        switch OrbCase

            case 'Lyapunov'

                % Run single shooting algorithm for Lyapunov and other planar
                % orbits symmetric about the x-axis
                [X,T_half,DF_fin,sshoot] = SingShoot_Lyap_PArc(X0,T0_half,sshoot);

            case 'Halo'
                
                % Run single shooting algorithm for Lyapunov and other planar
                % orbits symmetric about the x-axis
                [X,T_half,DF_fin,sshoot] = SingShoot_Halo_PArc(X0,T0_half,sshoot);
                
            case 'Vertical'
                
        end
        
    case 'NatParam'

        switch OrbCase

            case 'Lyapunov'

                % Run single shooting algorithm for Lyapunov and other planar
                % orbits symmetric about the x-axis
                [X,T_half,DF_fin,sshoot] = SingShoot_Lyap_NatParam(X0,T0_half,sshoot);

            case 'Halo'

            case 'Vertical'
                
        end
        
end