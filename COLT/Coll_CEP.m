function [Z,x_bnd,t_var,t_bnd,C,colt] = Coll_CEP(Z,t_bnd,t_var,header,colt)
% function [Z,x_bnd,t_var,t_bnd,C,colt] = Coll_CEP(Z,t_bnd,t_var,header,colt)
% 
% Given an initial design variable vector and a set of boundary node times 
% this function solves the collocation problem and implements CEP mesh 
% refinement.
%
% INPUTS:
%    Z         initial design variable vector ((n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  boundary nodes (n_seg x 1)
%    colt      structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    Z        final design variable vector (l*n_seg*(N+1)/2 x 1)
%    x_bnd    matrix of final boundary node states (l x 1 x n_seg+1)
%    t_var    non-normalized times at  final variable nodes (n_seg*(N+1)/2 x 1)
%    t_bnd    non-normalized times at  final boundary nodes (n_seg x 1)
%    C        matrix of polynomial coefficients (l x (N+1) x n_seg)
%
% Written by R. Pritchett, 10/20/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg_new = colt.n_seg; % initialize new segment value
n_seg_old = n_seg_new + 1; % initialize old segment value

% Initialize counters
colt.remi = 1; %initialize major while loop counter
colt.addi = 1; % initialize max error
colt.maxdiffe = 'N/A'; % no error distribution value for CEP mesh refinement
colt.maxe = 1; % initialize max error
colt.newti = 0; % initialize Newton's method iteration counter

% Print header for output table if switched on
if strcmp('On',header) == 1
    fprintf('\n                                    Sparsity                                                          \n')
    fprintf(' Deg  Segs  len(X)  len(F)  Non-0s    (%%%%)        max(|dX|)         max(|F|)       Iter   Rem   Add    max(|error|) \n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment Removal %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while n_seg_new < n_seg_old
    
    %Solve collocation problem
    header = 'Off'; % turn on iteration prinout
    [Z,x_bnd,tau_nodes,C,colt] = CollSolve_LT(Z,t_bnd,header,colt);
    
    % If max_iteration limit was reached exit Coll_CEP
    if colt.max_iter_chk 
        fprintf(2,'Warning: Coll_CEP stopped prematurely due to failure to converge collocation problem \n'); 
        return
    end

    %Convert column vector of design variables into 3D matrices
    [~,~,uis,sis,~,~] = Z23D(Z,colt);

    %Iterate counter
    colt.remi = colt.remi + 1;
    
    %Reset segment counter
    n_seg_old = n_seg_new;
    
    %Run mesh refinement function to remove segments
    [t_bnd_new,colt.maxe] = CEPrem(x_bnd,uis,t_bnd,colt);
        
    %Calculate new number of segments
    n_seg_new = size(t_bnd_new,1)-1;
    
    %Update segment number in colt structure
    colt.n_seg = n_seg_new;
    
    %Interpolate to obtain new variable nodes states
    [Z,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
        C,uis,sis,n_seg_old,n_seg_new,colt);
    
    %Reset mesh times
    t_bnd = t_bnd_new;
    t_var = t_var_new;
     
end

n_seg_old = n_seg_new-1; % adjust segment count in order to enter next loop
    
%Solve collocation problem once in between segment removal and addition loops
header = 'Off'; % turn on iteration prinout
[Z,x_bnd,tau_nodes,C,colt] = CollSolve_LT(Z,t_bnd,header,colt);

% If max_iteration limit was reached exit Coll_CEP
if colt.max_iter_chk 
    fprintf(2,'Warning: Coll_CEP stopped prematurely due to failure to converge collocation problem \n'); 
    return
end;

%Convert column vector of design variables into 3D matrices
[~,~,uis,sis,~,~] = Z23D(Z,colt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment Addition %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while n_seg_new > n_seg_old
    
    %Iterate counter
    colt.addi = colt.addi + 1;
    
    %Reset segment counter
    n_seg_old = n_seg_new;
    
    %Run mesh refinement function to remove segments
    [t_bnd_new,colt.maxe] = CEPadd(x_bnd,uis,t_bnd,colt);
        
    %Calculate new number of segments
    n_seg_new = size(t_bnd_new,1)-1;
    
    %Update segment number in colt structure
    colt.n_seg = n_seg_new;
    
    %Interpolate to obtain new variable nodes states
    [Z,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
        C,uis,sis,n_seg_old,n_seg_new,colt);
    
    %Reset mesh times
    t_bnd = t_bnd_new;
    t_var = t_var_new;
    
    %Solve collocation problem
    header = 'Off'; % turn on iteration prinout
    [Z,x_bnd,tau_nodes,C,colt] = CollSolve_LT(Z,t_bnd,header,colt);
    
    % If max_iteration limit was reached exit Coll_CEP
    if colt.max_iter_chk 
        fprintf(2,'Warning: Coll_CEP stopped prematurely due to failure to converge collocation problem \n'); 
        return
    end;
    
    %Convert column vector of design variables into 3D matrices
    [~,~,uis,sis,~,~] = Z23D(Z,colt);

end

%Final number of segments after addition and removal loops
colt.n_seg = n_seg_new; 

