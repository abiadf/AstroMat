function [t_nodes] = CollVarTimes(t_pts,N,n_seg,tau_nodes,mesh,pin)
% function [t_nodes] = CollVarTimes(t_pts,N,n_seg,tau_nodes,mesh,pin)
% 
% Given the necessary collocation parameters this function uses the 
% normalized times for nodes within a collocation segment to calculate the 
% non-normalized nondimensional times for the individual variable nodes 
% within each segment.  
%
% INPUTS:
%    t_pts      vector of nondimensional boundary node times 
%    N          degree of the interpolating polynomial    
%    n_seg      number of segments     
%    tau_nodes  vector of normalized variable and defect node times 
%               (N x 1) (same for all segments)      
%    mesh       indicates chosen node placement method, can be either LG 
%               or LGL (string data type)      
%    pin        indicates desired variation of LGL node placement method       
%
% OUTPUTS:
%    t_nodes    vector of non-normalized times for the variable nodes on
%               all segments (n_seg*(N+1)/2 x 1)
%
% Oringinally Written by: R. Pritchett, 07/24/2015
% Last Update: R. Pritchett, 09/24/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_nodes_var = tau_nodes(1:2:end); % select only variable (odd) nodes

switch mesh
    case 'LG'
        n_var = (N+1)/2; % number of variable nodes per segment
        t_nodes = zeros(n_seg*(N+1)/2,1); % preallocate
        %Convert from normalized to non-normalized time
        for ii = 1:n_seg
            t_nodes(n_var*ii - (n_var-1):n_var*ii) = ...
                0.5*(tau_nodes_var*(t_pts(ii+1) - t_pts(ii)) + t_pts(ii+1) + t_pts(ii));
        end
        
    case 'LGL'
        
        if pin
            n_var = (N+1)/2; % number of variable nodes per segment
            t_nodes = zeros(n_seg*(N+1)/2,1);
            %Convert from normalized to non-normalized time
            for ii = 1:n_seg
                t_nodes(n_var*ii-(n_var-1):n_var*ii) = ...
                    0.5*(tau_nodes_var*(t_pts(ii+1) - t_pts(ii)) + t_pts(ii+1) + t_pts(ii));
            end
        else          
            n_def = (N-1)/2; % number of defect nodes per segment        
            t_nodes = zeros(n_seg*(N-1)/2+1,1);
            %Convert from normalized to non-normalized time
            for ii = 1:n_seg
                t_nodes(n_def*ii-(n_def-1):n_def*ii) = ...
                    0.5*(tau_nodes_var(1:end-1)*(t_pts(ii+1) - t_pts(ii)) + t_pts(ii+1) + t_pts(ii));
            end        
            t_nodes(end)=t_pts(end);
        end
end
