function [Z,x_bnd,t_var,t_bnd,C,colt] = DirectTrans_InEq(Z,t_bnd,t_var,colt)
% function [Z,x_bnd,t_var,t_bnd,C] = DirectTrans(Z,t_bnd,t_var,colt)
%
% This function executes direct transcription with collocation. User inputs
% determine whether the function implements collocation with optimization
% (direct transcription) or collocation only (numerical integration). User
% input also determines the type of mesh refinement that is implemented
% with choices including de Boor, CEP, or no mesh refinement.
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
%
% Originally Written by: R. Pritchett, 10/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
OptMeth = colt.OptMeth;
Mesh = colt.Mesh;
opt_max_fevals = colt.opt_max_fevals;
opt_max_iter = colt.opt_max_iter;

%-------------------------------------------------------------------------%
    % Select Optimization Method %
%-------------------------------------------------------------------------%
    
switch OptMeth % Select case based on desired optimization method
    
    case 'NoOpt' % No optimization performed
          
    %---------------------------------------------------------------------%
    % Select Mesh Refinement Method %
    %---------------------------------------------------------------------%
        
    switch Mesh % Select case based on desired mesh refinement method 
        
        case 'NoMesh' % No mesh refinement
            
            % Collocation no mesh refinement
            header = 'On'; % turn on iteration prinout
            newti = 0; % initialize Newton's method iteration counter
            [Z,x_bnd,~,C,~] = CollSolve_LT(Z,t_bnd,header,colt,newti);
   
        case 'deBoor' % de Boor mesh refinement method
            
            %Collocation with de Boor Mesh Refinement
            header = 'On'; % turn on iteration prinout
            [Z,x_bnd,t_var,t_bnd,C,~,~,colt] = Coll_deBoor(Z,t_bnd,t_var,header,colt);
            
        case 'CEP' % CEP mesh refinement method
            
            % Collocation with Control Explicit Propagation
            header = 'On'; % turn on iteration prinout
            [Z,x_bnd,t_var,t_bnd,C,~,~,colt] = Coll_CEP(Z,t_bnd,t_var,header,colt);  
             
    end
    
    case 'fmincon' % Optimization performed with Matlab's fmincon function
    
    %---------------------------------------------------------------------%
    % Select Mesh Refinement Method %
    %---------------------------------------------------------------------%
    
        switch Mesh
            
            case 'NoMesh' % optimization with no mesh refinement
                
                % Calculate sparsity pattern of the Jacobian and Hessian matrices
                [collmat] = OptSetup(t_bnd,colt);
                
                % Define optimization method options
                options_opt = optimoptions('fmincon','Algorithm','interior-point',...
                    'CheckGradients',false,'Display','iter-detailed',...
                    'MaxFunctionEvaluations',opt_max_fevals,'MaxIterations',opt_max_iter,...
                    'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
                 
                % Run fmincon algorithm
                [Z,fval,exitflag,output,lambda,grad,Hoff] = ...
                    fmincon(@(Z) DrctTrans_Obj_MaxMf(Z,t_bnd,colt,collmat),Z,...
                    [],[],[],[],colt.lb_Z,colt.ub_Z,...
                    @(Z) DrctTrans_Con(Z,t_bnd,colt,collmat),options_opt);
                
                % Process fmincon output
                newti = 0; % initialize Newton's method iteration counter
                [x_bnd,~,C,~,colt] = post_fmincon(Z,t_bnd,output,newti,colt);

            case 'deBoor' % optimzation with de Boor mesh refinement
                
                [Z,x_bnd,t_var,t_bnd,C,~,~] = Opt_Coll_deBoor(Z,t_bnd,t_var,colt);
        
            case 'CEP'
            
            % To be coded at a later date
        
        end
        
    case 'IPOPT' % Optimization 
        
    % To be coded at a later date
   
end
        
        
