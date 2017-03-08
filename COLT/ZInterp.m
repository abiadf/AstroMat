function [Z,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,C,uis,sis,n_seg,n_seg_new,colt)
% function [Z,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,C,uis,sis,n_seg,n_seg_new,colt)
% 
% This function calculates the design variable vector for a new mesh given 
% a set of boundary nodes times for both the old and new mesh. It was 
% written for use in a collocation mesh refinement algorithm. 
%
% INPUTS:
%    t_bnd        non-normalized times at boundary nodes for original mesh (n_seg x 1)
%    t_bnd_new    non-normalized times at boundary nodes for new mesh (n_seg x 1)
%    t_var        non-normalized times at variable nodes for old mesh (n_seg_new*(N+1)/2 x 1)
%    tau_nodes    vector of normalized variable and defect node times (N x 1) (same for all segments)
%    C            matrix of polynomial coefficients (l x (N+1) x n_seg)
%    uis          3D matrix of control variable values (n_cntrl x (N+1)/2 x n_seg)
%    sis          3D matrix of slack variable values (n_slack x (N+1)/2 x n_seg)
%    n_seg        number of segments in old mesh
%    n_seg_new    number of segments in new mesh
%    colt         structure containing collocation and optimization parameters
%
% OUTPUTS:
%    X            design variable vector (l*n_seg_new*(N+1)/2 x 1)
%    t_var_new    non-normalized times at variable nodes for new mesh (n_seg_new*(N+1)/2 x 1)
%
% Written by R. Pritchett, 7/27/15
% Last Update: R. Pritchett, 09/29/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
NodeSpace = colt.NodeSpace;
N = colt.N;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
l = n_cntrl+n_slack+n_state*(N+1)/2;

%Calculate the non-normalized times of all variable nodes
t_var_new = CollVarTimes(t_bnd_new,N,n_seg_new,tau_nodes,NodeSpace,0);

% Reshape control and slack variable matrices 
uis_reshape = zeros(n_seg,n_cntrl);
sis_reshape = zeros(n_seg,n_slack);
for ii = 1:n_seg
   uis_reshape(ii,:) = uis(:,1,ii)';
   sis_reshape(ii,:) = sis(:,1,ii)';
end

%Use converged solution of the existing mesh to calculate X
cnt = 1; % initialize variable node counter variable
seg_cnt = 1; % initialize segment counter variable
x_col = zeros(n_state*(N+1)/2,1); % initialize column vector of state variable values for a single segment
Z_mat = zeros(l,n_seg_new); % initialize matrix of design variable vector values
for ii = 1:length(t_var_new)

    %Interpolate states from polynomial
    t = t_var_new(ii,1);
    [x_new] = StateInterp(t,C,t_bnd,N);
    x_col(n_state*(cnt-1)+1:n_state*cnt) = x_new;

    % Actions to perform when the last variable node of the segment is reached
    if cnt == (N+1)/2
        
        % Check which t_bnd current t is between and choose control values accordingly
        tbnd_less = find(t_bnd<t); % calculates indices of t_var less than t
        if isempty(tbnd_less) == 1; tbnd_less = 1; end;
        u_new = uis_reshape(tbnd_less(end),:)'; % use index of last t_bnd less than t 
        s_new = sis_reshape(tbnd_less(end),:)'; % use index of last t_bnd less than t 
        
        % Fill column of Z_mat matrix corresponding to current segment
        Z_mat(:,seg_cnt) = [u_new;s_new;x_col];
        
        % Reset variable node counter variable
        cnt = 0;
        % Iterate segment counter variable
        seg_cnt = seg_cnt + 1;
        % Reset column vector of state variable values for a single segment
        x_col = zeros(n_state*(N+1)/2,1);
        
    end
        
    % Iterate variable node counter variable
    cnt = cnt+1; % variable to track which variable node within a given segment the current interpolation is for
    
end

% Reshape Z_mat into a single column vector
Z = reshape(Z_mat,[l*n_seg_new 1]);
