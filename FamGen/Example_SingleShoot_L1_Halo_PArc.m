% Example_SingleShoot_L1_Halo_PArc.m
%
% This script uses multiple shooting and pseudo-arclength continuation to 
% compute the family of L1 Halo orbits 
%
% Originally Written by: R. Pritchett, 02/25/2017
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
sshoot.newt_tol = 1e-12; % tolerance for Newton's Method
sshoot.atten_tol = 1e-4; % tolerance for attenuation factor
sshoot.atten = 1; % attenuation factor for Newton's method step size 1/atten
sshoot.DiffType = 'Forward'; % finite difference method for numerical jacobian
sshoot.n_state = 6; % number of state variables
sshoot.n_slack = 0; % number of slack variables

%%% Targeter Case Inputs %%%
sshoot.ContCase = 'NoCont'; % target single orbit, no continuation
% sshoot.ContCase = 'PArc'; % target orbit for use in pseudo-arclength arclength continuation scheme
% sshoot.ContCase = 'NatParam'; % target orbit for use in natural parameter continuation
% sshoot.OrbCase = 'Lyapunov'; % target Lyapunov or other planar orbits symmetric about the x-axis
sshoot.OrbCase = 'Halo'; % target halo orbits
% sshoot.OrbCase = 'Vertical'; % target vertical orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % IC for Northern L1 Halo, taken from Grebow 2010
% X0 = [.8273 0 .0944 0 .2095 0]; % [nondimensional] Northern L1 Halo; 
% T0 = 2.7838; % [nondimensional]
% T0_half = T0/2; % [nondimensional]

% IC for Northern L1 Halo (Begin Planar)
X0 = [0.82339 0 0.00496 0 0.12674 0]; % [nondimensional] Northern L1 Halo; 
T0 = 2.74316; % [nondimensional]
T0_half = T0/2; % [nondimensional]

% % IC for Northern L2 Halo, taken from Grebow 2010
% X0 = [1.1690 0 .0979 0 -.1946 0]; % [nondimensional] Northern L2 Halo;
% T0 = 3.3314; % [nondimensional] 
% T0_half = T0/2; % [nondimensional] 

% Define integration timespan
tspan = [0 T0];

% Integrate initial guess
[x_init,t_init] = GslInteg_CR3BP(X0,tspan,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal
view([6,12])
% view([0,0])
% view([90,0])

% Plot Moon - realistic image, function from Bonnie Prado
bodyplot('Textures\','Moon',r_M,(1-mu)*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Earth - realistic image, function from Bonnie Prado
% bodyplot('Textures\','Earth',r_E,-mu.*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot initial state and orbit
plot3(X0(1)*l_ch,X0(2)*l_ch,X0(3)*l_ch,'*k');
plot3(x_init(:,1)*l_ch,x_init(:,2)*l_ch,x_init(:,3)*l_ch,':r');

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
max_step = 0.001;

% Define Z_prev from converged solution
Z_prev = [X(1); X(3); X(5); T_half];

% Calculate initial guess for next step of pseudo-arclenth continuation method
Zpls = Z_prev + delt_s.*delt_Z_prev;    

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
    X0_pls = [Zpls(1) X0(2) Zpls(2) X0(4) Zpls(3) X0(6)];
    T0_half_pls = Zpls(4);
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
    Z_prev = [X(1); X(3); X(5); T_half];
    
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
n_int = 1;

% Initialize iteration counter
iter = 0;

% Preallocate storage matrices
n_str = floor(n_step/n_int);
STM_store = zeros(6,6,n_str);
r_peri_store = zeros(n_str,1);
eig_store = zeros(n_str,6);
x0_store = zeros(n_str,6);
tf_store = zeros(n_str,1);
evec_store = zeros(6,6,n_str);

% Calculate Jacobi constant of each orbit
x0_nrho = reshape(cell2mat(X_store),[6 n_step])';
[jac_nrho] = Jacobi_Calc(x0_nrho(1:n_int:n_step,:),mu);

% Define a colormap
cmap = colormap(jet(n_str));

% Sort the Jacobi constant for colormap (so colors depend on the Jacobi constant of a particular curve)
[jac_sorted,jac_sorted_ind] = sort(jac_nrho);

for ii = 1:n_int:n_step
    
    iter = iter + 1;
    
    % Extract design variable vector for current step
    X0_i = X_store{ii};
    T_half_i = Thalf_store{ii};
    
    % Integrate converged result guess
    [x_fin,t_init] = GslInteg_STM([X0_i phi0],[0 2*T_half_i],mu);
    
    % Plot converged multiple shooting result
%     plot3(x_fin(1,1),x_fin(1,2),x_fin(1,3),'.r');
%     plot3(x_fin(end,1),x_fin(end,2),x_fin(end,3),'.g');
    plot3(x_fin(:,1).*l_ch,x_fin(:,2).*l_ch,x_fin(:,3).*l_ch,'b');
%     cmap_ind = find(iter == jac_sorted_ind);
%     plot3(x_fin(:,1).*l_ch,x_fin(:,2).*l_ch,x_fin(:,3).*l_ch,'color',cmap(cmap_ind,:));
    
    % Assemble and store the STM
    STM_store(:,:,iter) = reshape(x_fin(end,7:42),[6 6])';
    
    % Store orbit initial conditions and period
    x0_store(iter,:) = X0_i;
    tf_store(iter,:) = 2*T_half_i;
    
    % Propagate to periapse in order to calculate perilune radii
    tspan = [0 2*T_half_i]; 
    ic = X0_i;
    event_id = 2; % row vector that identifies the event type
    event_state = 0; % row vector defining desired event states
    dir = 0; % row vector indicating desired event crossing directions (+1 - positive crossing, -1 - negative crossing, 0 - all crossings).
    stop = 1; % row vector indicating whether integration should terminate at the events
    [x,~,~,~,~] = GslInteg_CR3BP_events(ic,tspan,mu,event_id,event_state,dir,stop);
    r_peri_store(iter,:) = norm(x(end,1:3) - [1-mu 0 0])*l_ch; % store perilune radius
    
    % Calculate and store eigenvalues and eigenvectors of monodromy matrix
    [evec_i,eig_i] = eig(STM_store(:,:,iter));
    eig_store(iter,:) = diag(eig_i);
    evec_store(:,:,iter) = evec_i;
    
end

% Label plot
xlabel('X [km]')
ylabel('Y [km]')
ylabel('Z [km]')
% caxis([min(jac_nrho) max(jac_nrho)]);
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
plot(r_peri_store,stable_ind_mag(:,1),'.b','MarkerSize',10);
plot(r_peri_store,stable_ubnd,':k','LineWidth',1);
plot(r_peri_store,stable_lbnd,':k','LineWidth',1);
% plot(Jac,stable_ind_mag(:,1),'.b','MarkerSize',10);
% plot(Jac,stable_ubnd,':k','LineWidth',1);
% plot(Jac,stable_lbnd,':k','LineWidth',1);
% xlabel('Jacobi Constant [nd]')
xlabel('Perilune Radius [km]')
ylabel('|\nu_1|')
title('L_1 Halo - Stability Index, \nu, as a Function of JC')
hold off

% Index 2
subplot(3,1,2)
hold on
grid on
plot(r_peri_store,stable_ind_mag(:,2),'.b','MarkerSize',10);
plot(r_peri_store,stable_ubnd,':k','LineWidth',1);
plot(r_peri_store,stable_lbnd,':k','LineWidth',1);
% plot(Jac,stable_ind_mag(:,2),'.b','MarkerSize',10);
% plot(Jac,stable_ubnd,':k','LineWidth',1);
% plot(Jac,stable_lbnd,':k','LineWidth',1);
% xlabel('Jacobi Constant [nd]')
xlabel('Perilune Radius [km]')
ylabel('|\nu_2|')
hold off

% Index 3
subplot(3,1,3)
hold on
grid on
plot(r_peri_store,stable_ind_mag(:,3),'.b','MarkerSize',10);
plot(r_peri_store,stable_ubnd,':k','LineWidth',1);
plot(r_peri_store,stable_lbnd,':k','LineWidth',1);
% plot(Jac,stable_ind_mag(:,3),'.b','MarkerSize',10);
% plot(Jac,stable_ubnd,':k','LineWidth',1);
% plot(Jac,stable_lbnd,':k','LineWidth',1);
% xlabel('Jacobi Constant [nd]')
xlabel('Perilune Radius [km]')
ylabel('|\nu_3|')
hold off

% Plot Jacobi constant as a function of perilune radius
figure (3)
hold on
grid on
plot(r_peri_store,Jac,'.b','MarkerSize',10);
xlabel('Perilune Radius [km]')
ylabel('Jacobi Constant [nd]')
title('L_1 Halo - Jacobi Constant as a Function of Perilune Radius')
hold off

% Plot Jacobi constant as a function of orbtial period
figure (4)
hold on
grid on
plot(tf_store.*t_ch./86400,Jac,'.b','MarkerSize',10);
xlabel('Orbital Period [days]')
ylabel('Jacobi Constant [nd]')
title('L_1 Halo - Jacobi Constant as a Function of Orbital Period')
hold off

% Plot Jacobi constant as a function of orbtial index
figure (5)
hold on
grid on
plot(Jac,tf_store.*t_ch./86400,'.b','MarkerSize',10);
xlabel('Jacobi Constant [nd]')
ylabel('Orbital Period [days]')
title('L_1 Halo - Orbital Period as a Function of Jacobi Constant')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select NRHO Orbits  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% NRHO Perilune Radii of Stability Changes (from Emily Zimovan) 
%-------------------------------------------------------------------------%

% Note from Emily: So, the way we are defining NRHOs is by 
% "bounded stability" meaning, although some of them are unstable, 
% they reside in some region where the stability remains "near by" the
% stable region (-1 <= nu <= 1).

% For the L1 NRHOs that includes the orbits with perilune radius 
% between 1135.86 km to 19062.5025 km.

% L1 NRHO Stability Change Altitudes:
% 934.34086363 km
% 1102.00521709 km
% 1135.8602697 km
% 6511.1936793 km
% 15682.4948721017 km
% 19062.5025232575 km

% Save IC for orbits with perilune radius between designated bounds
[nrho_ind] = find( 1136 < r_peri_store & r_peri_store < 19062);
x0_L1_nrho = x0_store(nrho_ind,:);
tf_L1_nrho = tf_store(nrho_ind);
n_L1_nrho = length(nrho_ind);

% For the L2 NRHOs that includes the orbits with perilune radius 
% between 1832.634 km to 17390.648 km

% L2 NRHO Stability Change Altitudes:
% 1832.6341858 km
% 13417.288694 km
% 17390.64858277 km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Periodic Orbit Family as .mat File  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('EM_L1_Halo_500_orbits','x0_store','tf_store')
% save('EM_L1_Halo_2000_orbits','x0_store','tf_store')
% save('EM_L1_NRHO_orbits','x0_L1_nrho','tf_L1_nrho','n_L1_nrho')

