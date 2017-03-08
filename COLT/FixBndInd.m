function [x0_des_ind,xf_des_ind] = FixBndInd(x0_des_ch,xf_des_ch)
% function [x0_des_ind,xf_des_ind] = FixBndInd(x0_des_ch,xf_des_ch)
% 
% This function calculates the constraint vector for collocation problems 
% involving low-thrust. 
%
% INPUTS:
%    x0_des_ch    cell array of character strings that designate whether
%                 each component of the initial boundary node is fixed or free
%    xf_des_ch    cell array of character strings that designate whether
%                 each component of the final boundary node is fixed or free
%
% OUTPUTS:
%    x0_des_ind    vector of indices of components of the initial boundary
%                  node that are fixed
%    xf_des_ind    vector of indices of components of the final boundary
%                  node that are fixed
%
% Written by R. Pritchett, 10/12/16
% Last Update: R. Pritchett, 10/12/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate number of components that are fixed or free
n_fix0 = length(x0_des_ch); 
n_fixf = length(xf_des_ch);

% Total number of components should be the same for both endpoints 
if n_fix0 ~= n_fixf
    msg = 'Number of components must be equal for initial and final endpoints';
    error(msg)
else
    n_fix = n_fix0;
end 

% Initialize iteration counters
cnt0 = 1;
cntf = 1;

% Initialize empty matrices
x0_des_ind = [];
xf_des_ind = [];

% Check each component of cell arrays to determine whether component is fixed or free
for ii = 1:n_fix
    
    % Define current component of x0 and xf cell arrays
    x0_des_chi = x0_des_ch{ii};
    xf_des_chi = xf_des_ch{ii};
    
    % If current component is fixed save its index value
    if strcmp(x0_des_chi,'fix') == 1
        x0_des_ind(cnt0) = ii;
        cnt0 = cnt0 + 1;
    end
    
    if strcmp(xf_des_chi,'fix') == 1
        xf_des_ind(cntf) = ii;
        cntf = cntf + 1;
    end
    
end
        