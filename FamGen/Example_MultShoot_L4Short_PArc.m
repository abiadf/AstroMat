% Example_MultShoot_L4Short_PArc.m
%
% This script uses multiple shooting and pseudo-arclength continuation to
% continue a family of L4 short period orbits (SPO). 
%
% Originally Written by: R. Pritchett, 02/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc;

%% Set Constants %%

%Call script that calculates useful CR3BP constants
CR3BPConst_EM

% Define collocation and optimization parameter structure
mult = struct;

% Save CR3BP constant parameters in structure 
mult.mu = mu; % gravitational 
mult.l_ch = l_ch; % characteristic length
mult.t_ch = t_ch; % characteristic time
mult.L = L; % libration point position coordinates 
mult.phi0 = phi0; % initial phi matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Numerical Integration Options %%%
RelativeTol = 1e-12; % relative integration tolerance
AbsoluteTol = 1e-12; % absolute integration tolerance
mult.options = odeset('RelTol',RelativeTol,'AbsTol',AbsoluteTol);
mult.atten_tol = 1e-4; % tolerance for attenuation factor
mult.atten = 1; % attenuation factor for Newton's method step size 1/atten

%%% Shooting Inputs %%%
mult.newt_tol = 1e-12; % tolerance for Newton's Method
mult.DiffType = 'Central'; % finite difference method for numerical jacobian
mult.n = 40; % number of patch points to use
mult.n_state = 6; % number of state variables
mult.n_slack = 0; % number of slack variables

%%% Boundary Constraints %%%

% Specify Boundary Constraint Case
% mult.BndCase = 'FixEndPtOnly';
mult.BndCase = 'Periodicity';

% Periodicity constraints
mult.period = [1 2 3 4 6]; % define which five states to include in periodicity constraint
mult.n_period = length(mult.period); % calculate number of periodicity constraints

% Initial Hyperplane Constraint
mult.i_hypr0 = [2]; % defines which component of the initial state is constrained to a hyperplane 
mult.x_hypr0 = [L(4,2)]; % defines the value(s) of the component(s) that define the hyperplane 
mult.n_hypr0 = length(mult.i_hypr0); % calculate number of hyperplane constraints

% Final Hyperplane Constraint
mult.i_hyprf = [2 4 6]; % defines which component of the final state is constrained to a hyperplane 
mult.x_hyprf = [0 0 0]; % defines the value(s) of the component(s) that define the hyperplane 
mult.n_hyprf = length(mult.i_hyprf); % calculate number of hyperplane constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define initial conditions
T0 = 6.4; % [nondimensional] initial guess at orbit period

% Discretize initial guess and assemble design variable vector
[Z0,x_ppt,t_ppt,x_plot] = L4Short_Discretize(T0,mult);

% Calculate length of design variable vector
mult.Zl = length(Z0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal

% Plot Moon - realistic image, function from Bonnie Prado
bodyplot('Textures\','Moon',r_M,(1-mu)*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Earth - realistic image, function from Bonnie Prado
bodyplot('Textures\','Earth',r_E,-mu*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot L_4 and Initial Guess
% plot(L(4,1),L(4,2),'*k'); 
plot(x_plot(1,1).*l_ch,x_plot(1,2).*l_ch,'*r');
plot(x_plot(:,1).*l_ch,x_plot(:,2).*l_ch,'.r');

% Label plot
xlabel('X [km]')
ylabel('Y [km]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converge Initial Solution and Calculate Null Space Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n             \n')
fprintf('Initial Iteration')

% Run multiple shooting algorithm
[Z,DF_fin,mult] = MultShoot(Z0,mult);

% Save or load collocation results
% save('MSLyap_n10','Z','x_bnd','t_var','t_bnd','C')
% load MSLyap_n10.mat

% Convert DF_fin from sparse to full matrix
DF_fin_full = full(DF_fin);

% Calculate null vector from collocation results
delt_Z_prev = null(DF_fin_full);

% Define pseudo-arclength continuation scaling parameter
delt_s = 0.05;

% Define max allowable step size
max_step = 0.5;

% Label Z from converged solution as Zprev, since new Z will be Zpls
Z_prev = Z;

% Calculate initial guess for next step of pseudo-arclenth continuation method
Zpls = Z_prev + delt_s.*delt_Z_prev;    

% Store pseudo-arclength parameters in a structure
parc = struct;
parc.Z_prev = Z_prev;
parc.deltZ_prev = delt_Z_prev;
parc.delt_s = delt_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Pseudo-Arclength Constraints and Step Along Tangent to the Solution Space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define number of pseudo-arclength continuation steps to take
n_step = 600;

% Preallocate cell arrays for storing results
Z_store = cell(1,n_step);
delt_s_store = zeros(n_step,1);

% Initialize null_chk structure
null_chk = struct;
null_chk.flip_cnt = 0; % Initialize sign flip counter 
null_chk.flip_iter = []; % Initialize sign flip iteration counter 
null_chk.bif_cnt = 0; % Initialize bifurcation counter

for ii = 1:n_step
    
    fprintf('\n             \n')
    fprintf('Iter:')
    fprintf(' %1.0f',ii)      
    
    % Current delt_Z_prev is now delt_Z_prev_old
    delt_Z_prev_old = delt_Z_prev;
 
    % Run multiple shooting algorithm
    [Z,DF_fin,mult] = MultShoot_PArc(Zpls,mult,parc);

    % Save or load collocation results
    % save('MSLyap_n10','Z','x_bnd','t_var','t_bnd','C')
    % load MSLyap_n10.mat

    % Store Results
    Z_store{ii} = Z;
    
    % Convert DF_fin from sparse to full matrix
    DF_fin_full = full(DF_fin);

    % Calculate null vector from collocation results
    delt_Z_prev = null(DF_fin_full);
    
    % Check sign and dimension of null space
    null_chk.Zpls = Zpls;
    null_chk.ps_ind = ii;
    [delt_Z_prev,null_chk] = NullCheck(delt_Z_prev,delt_Z_prev_old,null_chk);
        
    % Label Z from converged solution as Zprev, since new Z will be Zpls
    Z_prev = Z;
    
    % Adjust step-size based on number of iterations required to converge
%     if mult.newti > 5 % # of iter greater than 5 divide step size by 2
%         delt_s = delt_s/2;
%     elseif mult.newti <= 5 && 2*delt_s < max_step  % # of iter less than 5 and max step size is not reached multiply step size by 2
%         delt_s = 2.*delt_s;
%     elseif mult.newti <= 5 && 2*delt_s > max_step  % # of iter less than 5 and max step size is reached set step size equal to max step size
%         delt_s = max_step;
%     end
    
    % Store step size value
    delt_s_store(ii,1) = delt_s;
    
    % Calculate initial guess for next step of pseudo-arclenth continuation method
    Zpls = Z_prev + delt_s.*delt_Z_prev;        
    
    % Store pseudo-arclength parameters in a structure
    parc.Z_prev = Z_prev;
    parc.deltZ_prev = delt_Z_prev;
    parc.delt_s = delt_s;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Final Result Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define interval between plotting orbits, i.e. n_int = 5 means plot every fifth orbit computed
n_int = 1;

% Initialize iteration counter
iter = 0;

% Preallocate storage matrices
n_str = floor(n_step/n_int);
x0_store = zeros(n_str,6);
tf_store = zeros(n_str,1);
STM_store = zeros(6,6,n_str);
eig_store = zeros(n_str,6);
evec_store = zeros(6,6,n_str);

% Calculate Jacobi constant of each orbit
Z_mat = cell2mat(Z_store);
x0_spo = Z_mat(1:6,:)';
[jac_spo] = Jacobi_Calc(x0_spo(1:n_int:n_step,:),mu);

% Define a colormap
cmap = colormap(jet(n_str));

% Sort the Jacobi constant for colormap (so colors depend on the Jacobi constant of a particular curve)
[jac_sorted,jac_sorted_ind] = sort(jac_spo);

for ii = 1:n_int:n_step
    
    % Iterate iteration counter
    iter = iter + 1;
    
    % Extract design variable vector for current step
    Z_i = Z_store{ii};

    % Convert column vector of design variables into matrix of state variables
    x_ppt_i = reshape(Z_i(1:end-1),[mult.n_state mult.n]);

    % Calculate patch point times using total time
    T_i = Z_i(end);
%     t_ppt_i = linspace(0,T_i,mult.n);
    
%     % Convert patch points into continuous trajectory
%     [x_traj_i,t_traj_i] = ppt2Traj(x_ppt_i,t_ppt_i,mult);
 
    % Define orbit initial conditions for one revolution
    x0_rev = x_ppt_i(:,1)';
    tspan_rev = [0 T_i];
    
    % Propagate
    [x_rev,t_rev] = GslInteg_STM([x0_rev phi0],tspan_rev,mu);
    
    % Store orbit initial conditions and period
    x0_store(iter,:) = x0_rev;
    tf_store(iter,:) = T_i;
    
    % Assemble and store the STM
    STM_store(:,:,iter) = reshape(x_rev(end,7:42),[6 6])';
    
    % Calculate and store eigenvalues and eigenvectors of monodromy matrix
    [evec_i,eig_i] = eig(STM_store(:,:,iter));
    eig_store(iter,:) = diag(eig_i);
    evec_store(:,:,iter) = evec_i;
    
    % Plot converged multiple shooting result
%     plot(x_traj_i(:,1).*l_ch,x_traj_i(:,2).*l_ch,'b');
    plot(x_rev(:,1).*l_ch,x_rev(:,2).*l_ch,'b');
%     cmap_ind = find(iter == jac_sorted_ind);
%     plot(x_rev(:,1).*l_ch,x_rev(:,2).*l_ch,'color',cmap(cmap_ind,:));


    
end

% Label plot
xlabel('X [km]')
ylabel('Y [km]')
% caxis([min(jac_spo) max(jac_spo)]);
% h = colorbar;
% ylabel(h,'Jacobi Constant [nd]')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort Eigenvalues and Plot Stability Indices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort eigenvalues
[eig_sorted,evec_sorted] = SortEig(eig_store,evec_store);

% Calculate stability indices for each orbit; sum pairs and multiply by 0.5
stable_ind = 0.5.*(eig_sorted(:,1:2:6) + eig_sorted(:,2:2:6));
stable_ind_mag = sqrt(real(stable_ind).^2 + imag(stable_ind).^2);

% Calculate Jacobi constant for each orbit
[Jac] = Jacobi_Calc(x0_store,mu);

% Create array of orbit index numbers
orb_ind = [1:n_str]';

% Create upper and lower stability index bounds for plotting
stable_ubnd = ones(1,n_str);
stable_lbnd = -ones(1,n_str);

% Plot Stability indices 
figure (2)

% Index 1
subplot(3,1,1)
hold on
grid on
plot(Jac,stable_ind_mag(:,1),'.b','MarkerSize',10);
plot(Jac,stable_ubnd,':k','LineWidth',1);
plot(Jac,stable_lbnd,':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Orbital Period [km]')
ylabel('|\nu_1|')
title('L_4 SPO - Stability Index, \nu, as a Function of JC')
hold off

% Index 2
subplot(3,1,2)
hold on
grid on
plot(Jac,stable_ind_mag(:,2),'.b','MarkerSize',10);
plot(Jac,stable_ubnd,':k','LineWidth',1);
plot(Jac,stable_lbnd,':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Perilune Radius [km]')
ylabel('|\nu_2|')
hold off

% Index 3
subplot(3,1,3)
hold on
grid on
plot(Jac,stable_ind_mag(:,3),'.b','MarkerSize',10);
plot(Jac,stable_ubnd,':k','LineWidth',1);
plot(Jac,stable_lbnd,':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Perilune Radius [km]')
ylabel('|\nu_3|')
hold off

% % Plot Jacobi constant as a function of perilune radius
% figure (3)
% hold on
% grid on
% plot(r_peri_store,Jac,'.b','MarkerSize',10);
% xlabel('Perilune Radius [km]')
% ylabel('Jacobi Constant [nd]')
% title('Jacobi Constant as a Function of Perilune Radius')
% hold off

% Plot Jacobi constant as a function of orbtial period
figure (3)
hold on
grid on
plot(Jac,tf_store.*t_ch./86400,'.b','MarkerSize',10);
xlabel('Jacobi Constant [nd]')
ylabel('Orbital Period [days]')
title('L_4 SPO - Jacobi Constant as a Function of Orbital Period')
hold off

% Plot Jacobi constant as a function of orbital index
figure (4)
hold on
grid on
plot(orb_ind,Jac,'.b','MarkerSize',10);
xlabel('Orbital Index')
ylabel('Jacobi Constant [nd]')
title('L_4 SPO - Jacobi Constant as a Function of Orbital Index')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Periodic Orbit Family as .mat File  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save number of orbits computed
n_spo = length(tf_store);

% Save data
% save('EM_L4_SPO_600_orbits','x0_store','tf_store','n_spo')
    





