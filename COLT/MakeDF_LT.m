function [DF] = MakeDF_LT(Z,t_seg,t_seg_d,collmat,colt)
% function [DF] = MakeDF_LT(Z,t_seg,t_seg_d,collmat,colt)
%
% This function calculates the Jacobian matrix for collocation problems 
% involving low-thrust when coast parameters are included. The Jacobian 
% matrix in collocation problems is highly sparse, thus this function 
% calculates only the partials that are not guaranteed to be nonzero 
% values. However, some value in the output DF vector may still be zero.
%
% INPUTS:
%    Z          design variable vector (n_coast+(n_state+n_cntrl+n_slack)*n_seg*(N+1)/2 x 1)
%    t_seg      non-normalized times for each segment (l x (N+1)/2 x m)
%    t_seg_d    non-normalized times for each segment (l x (N-1)/2 x m)
%    collmat    structure containing constant collocation matrices
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    DF         nonzeros of the Jacobian matrix in column vector form
%
% Written by R. Pritchett, 6/07/16
% Last Update: R. Pritchett, 09/29/2016
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
Tmax = colt.Tmax;
x0_des_ind = colt.x0_des_ind;
xf_des_ind = colt.xf_des_ind;

% Calculate number of fixed components in initial and final endpoints
n_fix0 = length(x0_des_ind);
n_fixf = length(xf_des_ind);

% Convert column vector of design variables to 3D matrix
[zis,xis,uis,sis,coast_times,l] = Z23D(Z,colt);

%Set step size for complex step differentiation
stpsz = 1e-100; 

%Preallocate matrices of partial derivatives for collocation constraints
Dx0 = zeros(n_state,l,n_seg);
Dxf = zeros(n_state,l,n_seg);
Ddelta = zeros(n_state*(N-1)/2,l,n_seg);

%Preallocate matrices of partial derivatives for boundary and slack constraints
Dbnd0 = zeros(n_fix0,l); % derivatives of initial fixed endpoint constraints
Dbndf = zeros(n_fixf,l); % derivatives of final fixed endpoint constraints
Duhat = zeros(1,3,n_seg); % derivatives of uhat magnitude constraint
Dslack = zeros(n_slack,2,n_seg); % derivatives of slack variable constraints

%Preallocate perturbation matrices
u_pert = zeros(n_cntrl,(N+1)/2,n_seg);
x_pert = zeros(n_state,(N+1)/2,n_seg);

%Initialize iteration counters
count = 0;
uhat_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partials wrt Control and Slack Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------%
% Loop Through Control Design Variable Nodes %
%-----------------------------------------------------------------%
for ii = 1:n_cntrl
    
    % Iterate iteration counter
    count = count + 1;
    
    % Simultaneously perturb 1 state of 1 variable node in all segments
    u_pert(ii,:,:)=1i*stpsz; % complex step

    % Add perturbation to variable node states
    uis_p=uis+u_pert;

    % Compute collocation constraints with perturbed uis
    [~,delta_p,x0_p,xf_p,~] = Con_Defect(xis,uis_p,t_seg,t_seg_d,collmat,colt);            

    % Compute collocation constraint partials
    Dx0(:,count,:)=imag(x0_p)/stpsz;
    Dxf(:,count,:)=-imag(xf_p)/stpsz;
    Ddelta(:,count,:)=reshape(imag(delta_p)/stpsz,[n_state*(N-1)/2 1 n_seg]);

    % Compute endpoint constraint partials
    Dbnd0_col = imag(x0_p(1:n_state,:,1))/stpsz; % partials for initial endpoint constraint
    Dbnd0(:,count) = Dbnd0_col(x0_des_ind); % keep only partials for fixed components
    Dbndf_col = imag(xf_p(1:n_state,:,end))/stpsz; % partials for final endpoint constraint
    Dbndf(:,count) = Dbndf_col(xf_des_ind); % keep only partials for fixed components

    % Compute uhat magnitude constraint partials; if a thrust pointing control variable is currently selected 
    if ii > n_cntrl-3 % Control variable corresponding to u1 of thrust pointing
        uhat_count = uhat_count + 1; % iterate control iteration counter
        Duhat(:,uhat_count,:) = uis(ii,1,:)./sqrt(sum(uis(n_cntrl-2:n_cntrl,1,:).^2,1));
    end
        
    % Reset control variable perturbation matrix
    u_pert=0*u_pert;   
            
end % end control variable loop

%-----------------------------------------------------------------%
% Loop Through Slack Design Variable Nodes %
%-----------------------------------------------------------------%
for ii = 1:n_slack
    
    % Iterate iteration counter
    count = count + 1;
                
    % Compute Slack Variable Constraint Partials wrt Beta_Tlow
    Dslack(1,1,:) = ones(1,1,n_seg);
    Dslack(1,2,:) = -2.*sis(1,1,:);
    
    % Compute Slack Variable Constraint Partials wrt Beta_Thi
    Dslack(2,1,:) = -ones(1,1,n_seg);
    Dslack(2,2,:) = -2.*sis(2,1,:);
    
end % end slack variable loop
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partials wrt State Variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%-----------------------------------------------------------------%
% Loop Through Variable Nodes %
%-----------------------------------------------------------------%
for jj = 1:(N+1)/2
    
    %-----------------------------------------------------------------%
    % Loop Through Design Variables at Each Variable Node %
    %-----------------------------------------------------------------%
    for ii = 1:n_state
        
        count = count + 1; % iterate iteration counter
        
        % Simultaneously perturb 1 state of 1 variable node in all segments
        x_pert(ii,jj,:)=1i*stpsz; % complex step 

        % Add perturbation to variable node states
        xis_p=xis+x_pert;

        % Compute collocation constraints with perturbed xis
        [~,delta_p,x0_p,xf_p,~] = Con_Defect(xis_p,uis,t_seg,t_seg_d,collmat,colt);            

        % Compute collocation constraint partials 
        Dx0(:,count,:)=imag(x0_p)/stpsz;
        Dxf(:,count,:)=-imag(xf_p)/stpsz;
        Ddelta(:,count,:)=reshape(imag(delta_p)/stpsz,[n_state*(N-1)/2 1 n_seg]);

        % Compute coast continuity constraint partials
        Dbnd0_col = imag(x0_p(1:n_state,:,1))/stpsz; % partials for initial endpoint constraint
        Dbnd0(:,count) = Dbnd0_col(x0_des_ind); % keep only partials for fixed components
        Dbndf_col = imag(xf_p(1:n_state,:,end))/stpsz; % partials for final endpoint constraint
        Dbndf(:,count) = Dbndf_col(xf_des_ind); % keep only partials for fixed components

        x_pert=0*x_pert; % reset perturbation matrix   
    
    end % end design variable loop
    
end % end variable node loop
                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reshape components of DF into a single column vector 
DF =   [  reshape(Dx0(:,:,2:end),[n_state*l*(n_seg-1) 1]);
          reshape(Dxf(:,:,1:end-1),[n_state*l*(n_seg-1) 1]);
          reshape(Ddelta,[n_state*l*(N-1)*n_seg/2 1]);
          reshape(Dbnd0, [n_fix0*l 1]);
          reshape(Dbndf, [n_fixf*l 1]); 
          reshape(Duhat,[3*n_seg 1]); 
          reshape(Dslack,[n_slack*2*n_seg 1]) ];
      