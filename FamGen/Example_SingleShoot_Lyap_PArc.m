% Example_SingleShoot_Lyap_PArc.m
%
% This script uses single shooting and pseudo-arclength continuation to
% compute the family of L2 Lyapunov orbits in the Earth-Moon system.
%
% Originally Written by: R. Pritchett, 02/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc;

%% Set Constants %%

%Call script that calculates useful CR3BP constants
CR3BPConst_EM

% Define collocation and optimization parameter structure
sshoot = struct;

% Save CR3BP constant parameters in structure 
sshoot.mu = mu; % gravitational 
sshoot.l_ch = l_ch; % characteristic length
sshoot.t_ch = t_ch; % characteristic time
sshoot.L = L; % libration point position coordinates 
sshoot.phi0 = phi0; % initial phi matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User Inputs %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Numerical Integration Options %%%
RelativeTol = 1e-12; % relative integration tolerance
AbsoluteTol = 1e-12; % absolute integration tolerance
sshoot.options = odeset('RelTol',RelativeTol,'AbsTol',AbsoluteTol);
sshoot.options_event = odeset('Events',@xcross,'RelTol',1e-12,'AbsTol',1e-12);

%%% Shooting Inputs %%%
sshoot.newt_tol = 1e-10; % tolerance for Newton's Method
sshoot.atten_tol = 1e-4; % tolerance for attenuation factor
sshoot.atten = 1; % attenuation factor for Newton's method step size 1/atten
sshoot.DiffType = 'Forward'; % finite difference method for numerical jacobian
sshoot.n_state = 6; % number of state variables
sshoot.n_slack = 0; % number of slack variables

%%% Targeter Case Inputs %%%
sshoot.ContCase = 'NoCont'; % target single orbit, no continuation
% sshoot.ContCase = 'PArc'; % target orbit for use in pseudo-arclength arclength continuation scheme
% sshoot.ContCase = 'NatParam'; % target orbit for use in natural parameter continuation
sshoot.OrbCase = 'Lyapunov'; % target Lyapunov or other planar orbits symmetric about the x-axis
% sshoot.OrbCase = 'Halo'; % target halo orbits
% sshoot.OrbCase = 'Vertical'; % target vertical orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set of initial conditions from Tom Pavlak
load EM_L2_Lyap_250_orbits

% Define initial conditions
orb_num = 50;
X0 = u_0_save(orb_num,:);
T0 = t_f_save(orb_num);
T0_half = T0/2;
tspan = [0 T0];

% Integrate initial guess
[x_init,t_init] = GslInteg_CR3BP(X0,tspan,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal

plot(X0(1),X0(2),'*k');
plot(x_init(:,1),x_init(:,2),':r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converge Initial Solution and Calculate Null Space Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n             \n')
fprintf('Initial Iteration')

% Run single shooting algorithm
[X,T_half,DF_fin,sshoot] = SingShoot(X0,T0_half,sshoot);

% Save or load single shooting results
% save('SSLyap','X','T_half')
% load SSLyap.mat

% Calculate null vector from single shooting results
delt_Z_prev = null(DF_fin);

% Define pseudo-arclength continuation scaling parameter
delt_s = 0.001;

% Define max allowable step size
max_step = 0.005;

% Define Z_prev from converged solution
Z_prev = [X(1); X(5); T_half];

% Calculate initial guess for next step of pseudo-arclenth continuation method
Zpls = Z_prev - delt_s.*delt_Z_prev;    

% Store pseudo-arclength parameters in a structure
sshoot.Z_prev = Z_prev;
sshoot.deltZ_prev = delt_Z_prev;
sshoot.delt_s = delt_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Pseudo-Arclength Constraints and Step Along Tangent to the Solution Space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change targeter case
sshoot.ContCase = 'PArc'; % target orbit for use in pseudo-arclength arclength continuation scheme

% Define number of pseudo-arclength continuation steps to take
n_step = 200;

% Preallocate cell arrays for storing results
X_store = cell(1,n_step);
Thalf_store = cell(1,n_step);
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
 
    % Run single shooting algorithm
    X0_pls = [Zpls(1) X0(2:4) Zpls(2) X0(6)];
    T0_half_pls = Zpls(3);
    [X,T_half,DF_fin,sshoot] = SingShoot(X0_pls,T0_half_pls,sshoot);

    % Save or load single shooting results
    % save('SSLyap','X','T_half')
    % load SSLyap.mat

    % Store Results
    X_store{ii} = X;
    Thalf_store{ii} = T_half;
    
    % Convert DF_fin from sparse to full matrix
    DF_fin_full = full(DF_fin);

    % Calculate null vector from collocation results
    delt_Z_prev = null(DF_fin_full);
    
    % Check sign and dimension of null space
    null_chk.Zpls = Zpls;
    null_chk.ps_ind = ii;
    [delt_Z_prev,null_chk] = NullCheck(delt_Z_prev,delt_Z_prev_old,null_chk);

    % Define Z_prev from converged solution
    Z_prev = [X(1); X(5); T_half];
    
    % Adjust step-size based on number of iterations required to converge
    if sshoot.newti > 5 % # of iter greater than 5 divide step size by 2
        delt_s = delt_s/2;
    elseif sshoot.newti <= 5 && 2*delt_s < max_step  % # of iter less than 5 and max step size is not reached multiply step size by 2
        delt_s = 2.*delt_s;
    elseif sshoot.newti <= 5 && 2*delt_s > max_step  % # of iter less than 5 and max step size is reached set step size equal to max step size
        delt_s = max_step;
    end
    
    % Store step size value
    delt_s_store(ii,1) = delt_s;
    
    % Calculate initial guess for next step of pseudo-arclenth continuation method
    Zpls = Z_prev + delt_s.*delt_Z_prev;    

    % Store pseudo-arclength parameters in a structure
    sshoot.Z_prev = Z_prev;
    sshoot.deltZ_prev = delt_Z_prev;
    sshoot.delt_s = delt_s;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Final Result Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define interval between plotting orbits, i.e. n_int = 5 means plot every fifth orbit computed
n_int = 5;

% Initialize iteration counter
iter = 0;

% Preallocate storage matrices
n_str = floor(n_step/n_int);
x0_store = zeros(n_str,6);
tf_store = zeros(n_str,1);
STM_store = zeros(6,6,n_str);
eig_store = zeros(n_str,6);
xcross_store = zeros(n_str,1);
evec_store = zeros(6,6,n_str);

for ii = 1:n_int:n_step
    
    iter = iter + 1;
    
    % Extract design variable vector for current step
    X0_i = X_store{ii};
    T_half_i = Thalf_store{ii};
    
    % Integrate converged result guess
    [x_fin,t_init] = GslInteg_STM([X0_i phi0],[0 2*T_half_i],mu);
    
    % Plot converged multiple shooting result
    plot(x_fin(1,1),x_fin(1,2),'.r');
    plot(x_fin(end,1),x_fin(end,2),'.g');
    plot(x_fin(:,1),x_fin(:,2),'b');
    
    % Store orbit initial conditions and period
    x0_store(iter,:) = X0_i;
    tf_store(iter,:) = 2*T_half_i;
    
    % Assemble and store the STM
    STM_store(:,:,iter) = reshape(x_fin(end,7:42),[6 6])';
    
    % Calculate and store eigenvalues and eigenvectors of monodromy matrix
    [evec_i,eig_i] = eig(STM_store(:,:,iter));
    eig_store(iter,:) = diag(eig_i);
    xcross_store(iter,:) = x_fin(1,1); % store x crossing value for use in sorting function
    evec_store(:,:,iter) = evec_i;
    
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sort Eigenvalues and Plot Stability Indices %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort eigenvalues
[eig_sorted,evec_sorted] = SortEig(eig_store,evec_store);

% Calculate stability indices for each orbit; sum pairs and multiply by 0.5
stable_ind = 0.5.*(eig_sorted(:,1:2:6) + eig_sorted(:,2:2:6));

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
plot(orb_ind,stable_ind(:,1),'.b','MarkerSize',10);
plot(orb_ind,stable_ubnd,':k','LineWidth',1);
plot(orb_ind,stable_lbnd,':k','LineWidth',1);

hold off

% Index 2
subplot(3,1,2)
hold on
grid on
plot(orb_ind,stable_ind(:,2),'.b','MarkerSize',10);
plot(orb_ind,stable_ubnd,':k','LineWidth',1);
plot(orb_ind,stable_lbnd,':k','LineWidth',1);
hold off

% Index 3
subplot(3,1,3)
hold on
grid on
plot(orb_ind,stable_ind(:,3),'.b','MarkerSize',10);
plot(orb_ind,stable_ubnd,':k','LineWidth',1);
plot(orb_ind,stable_lbnd,':k','LineWidth',1);
hold off

% Plot Jacobi constant
figure (3)
hold on
grid on
plot(orb_ind,Jac,'.b','MarkerSize',10);
hold off