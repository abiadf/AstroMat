function [tau_nodes,A,A_inv,B,D,W] = CollSetup(N,mesh)
% function [tau_nodes, A, A_inv, B, D, W] = CollSetup(N,mesh)
% 
% Given the degree of the interpolating polynomial and the type of
% node placement method this function calculates the constant matrices
% used in collocation as well as the normalized times of the
% variable and defect nodes.
%
% INPUTS:
%    N          degree of the interpolating polynomial
%    mesh       indicates chosen node placement method, can be either LG 
%               or LGL (string data type)      
%
% OUTPUTS:
%    tau_nodes  vector of normalized variable and defect node times (N x 1) 
%               (same for all segments)
%    A          A matrix (N+1 x N+1)
%    A_inv      inverse of the A matrix (N+1 x N+1)
%    B          B matrix (N+1 x 2+(N-1)/2)
%    D          D matrix (N+1 x (N-1)/2)
%    W          W matrix ((N-1)/2 x 1)
%
% Last Update: R. Pritchett, 12/29/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the number of variable nodes per segment
n_var = (N+1)/2;

%Compute the placement and weights of the variable nodes based on the
%chosen placement method
switch mesh
    case 'LG'
        [lg_nodes, w] = lgwt(N, -1, 1);
        tau_nodes = sort(lg_nodes);
    case 'LGL'
        [lgl_nodes, w] = lglnodes(N-1);
        tau_nodes = sort(lgl_nodes);
end
tau_nodes_odd = tau_nodes(1:2:end);
tau_nodes_even = tau_nodes(2:2:end);

%Calculate the A matrix and its inverse
A = [(ones(N+1,1) * tau_nodes_odd.') .^ ((0:N).'*ones(1,n_var)) ((0:N).'*ones(1,n_var)) .* ((ones(N+1,1) * tau_nodes_odd.') .^ ((-1:N-1).'*ones(1,n_var)))];
A(isnan(A)) = 0;
A_inv = inv(A);

%Calculate the B matrix
B = [(-ones(N+1,1)) .^ ((0:N).') (ones(N+1,1) * tau_nodes_even.') .^ ((0:N).'*ones(1,(N-1)/2)) ones(N+1,1)];

%Calculate the D matrix
D = ((0:N).'*ones(1,(N-1)/2)) .* ((ones(N+1,1) * tau_nodes_even.') .^ ((-1:N-1).'*ones(1,(N-1)/2)));
D(isnan(D)) = 0;

%Calculate the W matrix (using only defect node weights)
W = diag(w(2:2:N-1));
