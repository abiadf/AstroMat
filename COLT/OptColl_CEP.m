function [Z,x_bnd,t_var,t_bnd,C,colt] = OptColl_CEP(Z,t_bnd,t_var,colt)
% function [Z,x_bnd,t_var,t_bnd,C] = OptColl_CEP(Z,t_bnd,t_var,colt)
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
%    colt      structure containing collocation and optimization parameters      
%
% Written by R. Pritchett, 10/20/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg_new = colt.n_seg; % initialize new segment value
n_seg_old = n_seg_new + 1; % initialize old segment value
options_opt = colt.options_opt; % define optimization options

% Initialize counters
colt.remi = 1; %initialize major while loop counter
colt.addi = 1; % initialize max error
colt.maxdiffe = 'N/A'; % no error distribution value for CEP mesh refinement
colt.maxe = 1; % initialize max error
newti = 0; % initialize Newton's method iteration counter

% Initialize matrix of lower and upper bounds on design variables in single segment 
lb_i = [0;-Inf.*ones(3,1); -Inf.*ones(colt.n_state*(colt.N+1)/2,1,1)];
ub_i = [colt.Tmax;Inf.*ones(3,1); Inf.*ones(colt.n_state*(colt.N+1)/2,1,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment Removal %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while n_seg_new < n_seg_old
    
    % Calculate sparsity pattern of the Jacobian and Hessian matrices
    [collmat] = OptSetup(t_bnd,colt);

    % Run fmincon algorithm
    [Z,~,~,output,~,~,~] = ...
        fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
        [],[],[],[],colt.lb_Z,colt.ub_Z,...
        @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

    % Process fmincon output
    [x_bnd,tau_nodes,C,newti,colt] = PostOpt(Z,t_bnd,output,newti,colt);
    
    % Print output for one iteration of CEP mesh refinement
    Zl = length(Z); % calculate length of design variable vector
    Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
    error = output.constrviolation; % extract error from fmincon output
    fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Rem   Add    max(|error|) \n')
    fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e \n',N,n_seg,Zl,Fl,error,newti,colt.remi,colt.addi,colt.maxe)
    
    % Convert column vector of design variables into 3D matrices
    [~,~,uis,sis,~,l] = Z23D(Z,colt);

    % Iterate counter
    colt.remi = colt.remi + 1;
    
    % Reset segment counter
    n_seg_old = n_seg_new;
    
    % Run mesh refinement function to remove segments
    [t_bnd_new,colt.maxe] = CEPrem(x_bnd,uis,t_bnd,colt);
        
    % Calculate new number of segments
    n_seg_new = size(t_bnd_new,1)-1;
    
    % Update segment number in colt structure
    colt.n_seg = n_seg_new;
    
    % Interpolate to obtain new variable nodes states
    [Z,t_var_new] = ZInterp(t_bnd,t_bnd_new,t_var,tau_nodes,...
        C,uis,sis,n_seg_old,n_seg_new,colt);
    
    % Formulate lower and upper bounds on state variables
    lb_mat = repmat(lb_i,[1 1 colt.n_seg]);
    ub_mat = repmat(ub_i,[1 1 colt.n_seg]);
    colt.lb_Z = reshape(lb_mat,[l*colt.n_seg 1]);
    colt.ub_Z = reshape(ub_mat,[l*colt.n_seg 1]);
    
    % Reset mesh times
    t_bnd = t_bnd_new;
    t_var = t_var_new;
    
end

n_seg_old = n_seg_new-1; % adjust segment count in order to enter next loop
    
% Calculate sparsity pattern of the Jacobian and Hessian matrices
[collmat] = OptSetup(t_bnd,colt);

% Run fmincon algorithm
[Z,~,~,output,~,~,~] = ...
    fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
    [],[],[],[],colt.lb_Z,colt.ub_Z,...
    @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

% Process fmincon output
[x_bnd,tau_nodes,C,newti,colt] = PostOpt(Z,t_bnd,output,newti,colt);

% Print output for one iteration of CEP mesh refinement
Zl = length(Z); % calculate length of design variable vector
Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
error = output.constrviolation; % extract error from fmincon output
fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Rem   Add    max(|error|) \n')
fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e \n',N,n_seg,Zl,Fl,error,newti,colt.remi,colt.addi,colt.maxe)

%Convert column vector of design variables into 3D matrices
[~,~,uis,sis,~,l] = Z23D(Z,colt);

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
    
    % Formulate lower and upper bounds on state variables
    lb_mat = repmat(lb_i,[1 1 colt.n_seg]);
    ub_mat = repmat(ub_i,[1 1 colt.n_seg]);
    colt.lb_Z = reshape(lb_mat,[l*colt.n_seg 1]);
    colt.ub_Z = reshape(ub_mat,[l*colt.n_seg 1]);
    
    %Reset mesh times
    t_bnd = t_bnd_new;
    t_var = t_var_new;
    
    % Calculate sparsity pattern of the Jacobian and Hessian matrices
    [collmat] = OptSetup(t_bnd,colt);

    % Run fmincon algorithm
    [Z,~,~,output,~,~,~] = ...
        fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
        [],[],[],[],colt.lb_Z,colt.ub_Z,...
        @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);

    % Process fmincon output
    [x_bnd,tau_nodes,C,newti,colt] = PostOpt(Z,t_bnd,output,newti,colt);
    
    % Print output for one iteration of CEP mesh refinement
    Zl = length(Z); % calculate length of design variable vector
    Fl = length(lambda.eqnonlin); % calculate length of constraint vector 
    error = output.constrviolation; % extract error from fmincon output
    fprintf(' Deg  Segs  len(X)  len(F)     max(|F|)     Iter   Rem   Add    max(|error|) \n')
    fprintf(' %2.0f  %4.0f  %6.0f  %6.0f  %16.8e  %5.0f  %5.0f  %5.0f  %16.8e \n',N,n_seg,Zl,Fl,error,newti,colt.remi,colt.addi,colt.maxe)
    
    %Convert column vector of design variables into 3D matrices
    [~,~,uis,sis,~,l] = Z23D(Z,colt);

end

%Final number of segments after addition and removal loops
colt.n_seg = n_seg_new; 

