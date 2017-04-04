function [Z,x_bnd,t_var,t_bnd,C,colt] = OptColl_deBoor(Z,t_bnd,t_var,colt)
% function [Z,x_bnd,t_var,t_bnd,C,colt] = OptColl_deBoor(Z,t_bnd,t_var,colt)
% 
% Given an initial design variable vector and a set of boundary node times 
% this function solves the collocation problem and implements de Boor mesh 
% refinement.
%
% INPUTS:
%    Z         design variable vector (l*n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  boundary nodes (n_seg x 1)
%    t_var     non-normalized times at  variable nodes (n_seg*(N+1)/2 x 1)
%    colt      structure containing collocation and optimization parameters      
%
% OUTPUTS:
%    Z         design variable vector (l*n_seg*(N+1)/2 x 1)
%    x_bnd     matrix of final boundary node states (l x 1 x n_seg+1)
%    t_var     non-normalized times at  final variable nodes (n_seg*(N+1)/2 x 1)
%    t_bnd     non-normalized times at  final boundary nodes (n_seg x 1)
%    C         matrix of polynomial coefficients (l x (N+1) x n_seg)
%    colt      structure containing collocation and optimization parameters      
%
% Written by R. Pritchett, 09/24/16
% Last Update: R. Pritchett, 10/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
N = colt.N;
n_seg = colt.n_seg;
n_state = colt.n_state;
n_coast = colt.n_coast;
ctol = colt.ctol;
Dec = colt.Dec;
options_opt = colt.options_opt;

% Initialize counters
colt.maji = 1; %initialize major while loop counter
colt.maxe = 1; % initialize max error
colt.newti = 0; % initialize Newton's method iteration counter

% Initialize matrix of lower and upper bounds on design variables in single segment 
lb_i = [0;-Inf.*ones(3,1); -Inf.*ones(colt.n_state*(colt.N+1)/2,1,1)];
ub_i = [colt.Tmax;Inf.*ones(3,1); Inf.*ones(colt.n_state*(colt.N+1)/2,1,1)];

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
                
        % Calculate sparsity pattern of the Jacobian and Hessian matrices
        [collmat] = OptSetup(t_bnd,colt);

        % Run fmincon algorithm
        [Z,~,~,output,lambda,~,~] = ...
            fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
            [],[],[],[],colt.lb_Z,colt.ub_Z,...
            @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

        % Process fmincon output
        [x_bnd,tau_nodes,C,colt] = PostOpt(Z,t_bnd,output,colt);

        % Print statement for de Boor mesh refinement with fmincon
        Zl = length(Z); % calculate length of design variable vector
        Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
        error = output.constrviolation; % extract error from fmincon output
        fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Mesh   Node    max(|error|)     DErr \n')
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e  %5.0f \n',colt.N,colt.n_seg,Zl,Fl,error,colt.newti,colt.mini,colt.maji,colt.maxe,colt.maxdiffe)

        % Convert column vector of design variables into 3D matrices
        [~,~,uis,sis,~,~] = Z23D(Z,colt);

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
        
        % Add tau and alpha parameters to interpolated Z values
        Z_new = [Z(1:n_coast); Z_interp];
        
        % Formulate lower and upper bounds on state variables
        lb_mat = repmat(lb_i,[1 1 colt.n_seg]);
        ub_mat = repmat(ub_i,[1 1 colt.n_seg]);
        colt.lb_Z = reshape(lb_mat,[l*colt.n_seg 1]);
        colt.ub_Z = reshape(ub_mat,[l*colt.n_seg 1]);
               
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
        
        % Calculate sparsity pattern of the Jacobian and Hessian matrices
        [collmat] = OptSetup(t_bnd,colt);

        % Run fmincon algorithm
        [Z,~,~,output,lambda,~,~] = ...
            fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
            [],[],[],[],colt.lb_Z,colt.ub_Z,...
            @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

        % Process fmincon output
        [x_bnd,tau_nodes,C,colt] = PostOpt(Z,t_bnd,output,colt);
                
        % Print statement for de Boor mesh refinement with fmincon
        Zl = length(Z); % calculate length of design variable vector
        Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
        error = output.constrviolation; % extract error from fmincon output
        fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Mesh   Node    max(|error|)     DErr \n')
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e  %5.0f \n',colt.N,colt.n_seg,Zl,Fl,error,colt.newti,colt.mini,colt.maji,colt.maxe,colt.maxdiffe)

    else % if max error is above tolerance
        
        colt.maji = colt.maji + 1; %iterate max error loop iteration counter 
        
        %Update number of segments in mesh
        [n_seg_new,t_bnd_new] = UpdateMesh(t_bnd,I,xi,ctol,colt.maxe,N,n_seg);
        
        %Interpolate to obtain new variable node states
        [Z_interp,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
            C,uis,sis,n_seg,n_seg_new,colt);
               
        % Add tau and alpha parameters back to interpolated Z value 
        Z_new = [Z(1:n_coast); Z_interp];
        
        % Formulate lower and upper bounds on state variables
        lb_mat = repmat(lb_i,[1 1 colt.n_seg]);
        ub_mat = repmat(ub_i,[1 1 colt.n_seg]);
        colt.lb_Z = reshape(lb_mat,[l*colt.n_seg 1]);
        colt.ub_Z = reshape(ub_mat,[l*colt.n_seg 1]);
        
        % Reset mesh times and boundary node number
        t_bnd = t_bnd_new;
        t_var = t_var_new;
        Z = Z_new;
        n_seg = n_seg_new; 
        colt.n_seg = n_seg_new; % save new segment number in colt structure
        
        % Calculate sparsity pattern of the Jacobian and Hessian matrices
        [collmat] = OptSetup(t_bnd,colt);

        % Run fmincon algorithm
        [Z,~,~,output,lambda,~,~] = ...
            fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
            [],[],[],[],colt.lb_Z,colt.ub_Z,...
            @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

        % Process fmincon output
        [x_bnd,tau_nodes,C,colt] = PostOpt(Z,t_bnd,output,colt);
          
        % Print statement for de Boor mesh refinement with fmincon
        Zl = length(Z); % calculate length of design variable vector
        Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
        error = output.constrviolation; % extract error from fmincon output
        fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Mesh   Node    max(|error|)     DErr \n')
        fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e  %5.0f \n',colt.N,colt.n_seg,Zl,Fl,error,colt.newti,colt.mini,colt.maji,colt.maxe,colt.maxdiffe)

        % Convert column vector of design variables into 3D matrices
        [~,~,uis,sis,~,~] = Z23D(Z,colt);
        
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
