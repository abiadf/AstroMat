function [Z,x_bnd,t_var,t_bnd,C,n_seg,maxe,colt] = Coll_deBoor(Z,t_bnd,t_var,header,colt)
% function [Z,x_bnd,t_var,t_bnd,C,n_seg,maxe] = Coll_deBoor(Z,t_bnd,t_var,header,colt)
% 
% Given an initial design variable vector and a set of boundary node times 
% this function solves the collocation problem and implements de Boor mesh 
% refinement.
%
% INPUTS:
%    Z         design variable vector (l*n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  boundary nodes (n_seg x 1)
%    t_bnd     non-normalized times at  boundary nodes (n_seg x 1)
%    header    string that indicates whether to print output of each Newton's method iteration 
%    colt      structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    Z         design variable vector (l*n_seg*(N+1)/2 x 1)
%    x_bnd     matrix of final boundary node states (l x 1 x n_seg+1)
%    t_var     non-normalized times at  final variable nodes (n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  final boundary nodes (n_seg x 1)
%    C         matrix of polynomial coefficients (l x (N+1) x n_seg)
%    n_seg     final number of segments
%    maxe      max segment error
%
% Written by R. Pritchett, 09/24/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N;
n_seg = colt.n_seg;
n_cntrl = colt.n_cntrl;
n_state = colt.n_state;
n_slack = colt.n_slack;
n_coast = colt.n_coast;
ctol = colt.ctol;
Dec = colt.Dec;

% Initialize counters
colt.maji = 1; %initialize major while loop counter
colt.maxe = 1; % initialize max error
newti = 0; % initialize Newton's method iteration counter

% Print header for output table if switched on
if strcmp('On',header) == 1
    fprintf('\n                                    Sparsity                                                          \n')
    fprintf(' Deg  Segs  len(X)  len(F)  Non-0s    (%%%%)        max(|dX|)         max(|F|)       Iter   Mesh   Node    max(|error|)     DErr \n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collocation and de Boor Mesh Refinement %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Error Reduction Loop %
%-------------------------------------------------------------------------%

while colt.maxe > ctol % error reduction loop

    colt.mini = 1; % initialize distribution while loop counter

    %---------------------------------------------------------------------%
    % Error Distribution Loop %
    %---------------------------------------------------------------------%
    
    while colt.maxdiffe > 0 && colt.maxe > ctol % error distribution loop
                
        % Initial call to solve collocation problem
        header = 'Off'; % turn on iteration prinout
        [Z,x_bnd,tau_nodes,C,newti] = CollSolve_LT(Z,t_bnd,header,colt,newti);
        
        % Convert column vector of design variables into 3D matrices
        [zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);

        % Iterate distribution loop iteration counter
        colt.mini = colt.mini + 1;
        
        % Apply deBoor error equidistribution method
        [t_bnd_new,I,xi] = deBoor(t_bnd,C,N,n_state,n_seg);
        
        % Compute maxe and maxdiffe
        [colt.maxe,colt.maxdiffe] = ParseErr(t_bnd,xi,Dec,N,n_seg);
        
        % Interpolate to obtain new variable node states
        n_seg_new = n_seg; % number of segments not yet
        [Z_interp,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
            C,uis,sis,n_seg,n_seg_new,colt);
               
        %Add tau and alpha parameters to interpolated Z values
        Z_new = [Z(1:n_coast); Z_interp];
        
        % Old mesh and interpolated states are used if error distribution loop is exited 
        t_bnd_old = t_bnd;
        t_var_old = t_var;
        Z_old = Z;
        
        % New mesh and interpolated states are used if error distribution loop is repeated
        t_bnd = t_bnd_new;
        t_var = t_var_new;
        Z = Z_new;
    
    end
    
    %---------------------------------------------------------------------%
    % End - Error Distribution Loop %
    %---------------------------------------------------------------------%
    
    %%% Important Note: The calculated error values correspond to the  %%%
    %%% "old" mesh, be sure to use the "old" mesh when the mesh         %%%
    %%% distribution loop is exited                                    %%%
    
    % Old mesh and interpolated states are used 
    t_bnd = t_bnd_old;
    t_var = t_var_old;
    Z = Z_old;
    
    %---------------------------------------------------------------------%
    % Check Max Error %
    %---------------------------------------------------------------------%
    
    if colt.maxe < ctol % if max error is below tolerance
        
        %Solve collocation problem
        header = 'Off'; % turn on iteration prinout
        [Z,x_bnd,tau_nodes,C,newti] = CollSolve_LT(Z,t_bnd,header,colt,newti);

    else % if max error is above tolerance
        
        colt.maji = colt.maji + 1; %iterate max error loop iteration counter 
        
        %Update number of segments in mesh
        [n_seg_new,t_bnd_new] = UpdateMesh(t_bnd,I,xi,ctol,colt.maxe,N,n_seg);
        
        %Interpolate to obtain new variable node states
        [Z_interp,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
            C,uis,sis,n_seg,n_seg_new,colt);
               
        % Add coast parameters back to interpolated Z value 
        Z_new = [Z(1:n_coast); Z_interp];
        
        %Reset mesh times and boundary node number
        t_bnd = t_bnd_new;
        t_var = t_var_new;
        Z = Z_new;
        n_seg = n_seg_new; 
        colt.n_seg = n_seg_new; % save new segment number in colt structure
        
        % Solve collocation problem
        header = 'Off'; % turn on iteration prinout
        [Z,x_bnd,tau_nodes,C,newti] = CollSolve_LT(Z,t_bnd,header,colt,newti);

        % Convert column vector of design variables into 3D matrices
        [zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);
        
        %Apply deBoor error equidistribution method
        [~,I,xi] = deBoor(t_bnd,C,N,n_state,n_seg);
        
        %Compute maxe and maxdiffe
        [colt.maxe,colt.maxdiffe] = ParseErr(t_bnd,xi,Dec,N,n_seg);
        
        % Set old mesh values equal to current mesh in case error distribution loop is not rentered
        t_bnd_old = t_bnd;
        t_var_old = t_var;
        Z_old = Z;
        
    end
   
end     

% Define maxe variable for output
maxe = colt.maxe;
     