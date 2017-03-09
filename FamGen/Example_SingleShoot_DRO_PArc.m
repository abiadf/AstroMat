% Example_Lyapunov_EM.m
%
% This script uses collocation and direct optimization to compute a
% Lyapunov orbit 
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
sshoot.OrbCase = 'Lyapunov'; % target Lyapunov or other planar orbits symmetric about the x-axis
% sshoot.OrbCase = 'Halo'; % target halo orbits
% sshoot.OrbCase = 'Vertical'; % target vertical orbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% Initial Guess: Method 1
%-------------------------------------------------------------------------%

% % In this example an initial guess for a DRO is obtained from a 2BP 
% % approximation about the Moon.

% % Define orbit altitude above moon
% r_alt = 1000; % [km] altitude above Moon's surface
% r_dim = r_alt + r_M; % [km] distance from the center of the Moon
% 
% % Calculate orbital velocity
% v_dim = sqrt(mu_moon/r_dim);
% 
% % Assemble initial states
% X0_dim = [r_dim 0 0 0 v_dim 0];
% X0 = [(1-mu)+r_dim./l_ch 0 0 0 v_dim./v_ch 0];
% 
% % Calculate orbital period
% T0_dim = 2*pi*sqrt(r_dim^3/mu_moon);
% T0 = T0_dim./t_ch;
% T0_half = T0/2;

%-------------------------------------------------------------------------%
% Initial Guess: Method 2
%-------------------------------------------------------------------------%

% % Initial guess from Henon Paper
% Gamma = 3; % given
% xi = -0.25071;  % given
% T0_half = 0.35696;  % given
% x = 1-mu+xi*mu^(1/3); % eqn from Henon paper
% C = 3+Gamma*mu^(2/3); % eqn from Henon paper
% d = sqrt((x+mu)^2);
% r = sqrt((x-1+mu)^2);
% U = (1-mu)/d + mu/r;
% Ustr = U + 0.5*(x^2);
% y_dot = sqrt(2*Ustr-C);
% y_dot_alt = sqrt(3*xi^2+2/sqrt(xi^2)-Gamma);
% X0 = [x 0 0 0 y_dot 0];
% T0 = 2*T0_half;

%-------------------------------------------------------------------------%
% Initial Guess: Method 3
%-------------------------------------------------------------------------%

% Alternate initial guess, obtained by continuing to lower lunar altitude
% DRO's and saving IC. Used Henon guess to perform this continuation method
X0 = [0.982776 0 0 0 1.552646 0];
T0 = 0.020533;
T0_half = T0/2;

% Define timespan for integration
tspan = [0 T0];

% Integrate initial guess
[x_init,t_init] = GslInteg_CR3BP(X0,tspan,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal

% Plot Moon - realistic image, function from Bonnie Prado
bodyplot('Textures\','Moon',r_M,(1-mu).*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Earth - realistic image, function from Bonnie Prado
bodyplot('Textures\','Earth',r_E,-mu.*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot initial guess
plot3(x_init(1,1).*l_ch,x_init(1,2).*l_ch,x_init(1,3).*l_ch,'*r');
plot3(x_init(:,1).*l_ch,x_init(:,2).*l_ch,x_init(:,3).*l_ch,'.r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converge Initial Solution and Calculate Null Space Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run single shooting algorithm
[X,T_half,DF_fin,sshoot] = PeriodicTarget(X0,T0_half,sshoot);

% Save or load single shooting results
% save('SSLyap','X','T_half')
% load SSLyap.mat

% Calculate null vector from single shooting results
delt_Z_prev = null(DF_fin);

% Define pseudo-arclength continuation scaling parameter
delt_s = 0.005;

% Define max allowable step size
max_step = 0.02;

% Label Z from converged solution as Zprev, since new Z will be Zpls
Z_prev = [X(1); X(5); T_half];

% Calculate initial guess for next step of pseudo-arclenth continuation method
% delt_Z_prev = -1*delt_Z_prev;
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
n_step = 400;

% Preallocate cell arrays for storing results
X_store = cell(1,n_step);
Thalf_store = cell(1,n_step);
delt_s_store = zeros(n_step,1);
delt_Z_prev_store = zeros(3,n_step);
Zpls_store = zeros(3,n_step);

% Initialize null_chk structure
null_chk = struct;
null_chk.flip_cnt = 0; % Initialize sign flip counter 
null_chk.flip_iter = []; % Initialize sign flip iteration counter 
null_chk.bif_cnt = 0; % Initialize bifurcation counter

for ii = 1:n_step
    
    % Current delt_Z_prev is now delt_Z_prev_old
    delt_Z_prev_old = delt_Z_prev;
 
    % Run single shooting algorithm
    X0_pls = [Zpls(1) X0(2:4) Zpls(2) X0(6)];
    T0_half_pls = Zpls(3);
    [X,T_half,DF_fin,sshoot] = PeriodicTarget(X0_pls,T0_half_pls,sshoot);

    % Save or load single shooting results
    % save('SSLyap','X','T_half')
    % load SSLyap.mat

    % Store Results
    X_store{ii} = X;
    Thalf_store{ii} = T_half;
    
    % Calculate null vector from collocation results
    delt_Z_prev = null(DF_fin);
    
    % Check sign and dimension of null space
    null_chk.Zpls = Zpls;
    null_chk.ps_ind = ii;
%     [delt_Z_prev,null_chk] = NullCheck(delt_Z_prev,delt_Z_prev_old,null_chk);
    [delt_Z_prev,null_chk] = NullCheck_New(delt_Z_prev,delt_Z_prev_old,null_chk);
    
    % Label Z from converged solution as Zprev, since new Z will be Zpls
    Z_prev = [X(1); X(5); T_half];
    
    % Store step size value
    delt_s_store(ii,1) = delt_s;
    delt_Z_prev_store(:,ii) = delt_Z_prev;
    
    % Adjust step-size based on number of iterations required to converge
    if sshoot.newti > 10 % # of iter greater than 5 divide step size by 2
        delt_s = delt_s/2;
    elseif sshoot.newti <= 10 && 2*delt_s <= max_step  % # of iter less than 5 and max step size is not reached multiply step size by 2
        delt_s = 2.*delt_s;
    elseif sshoot.newti <= 10 && 2*delt_s > max_step  % # of iter less than 5 and max step size is reached set step size equal to max step size
        delt_s = max_step;
    end

    % Calculate initial guess for next step of pseudo-arclenth continuation method
    Zpls = Z_prev + delt_s.*delt_Z_prev;  
    
    % Store Z_pls at each iteration
    Zpls_store(:,ii) = Zpls; 

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
r_peri_store = zeros(n_str,1);
eig_store = zeros(n_str,6);
evec_store = zeros(6,6,n_str);

% Calculate Jacobi constant of each orbit
x0_dro = reshape(cell2mat(X_store),[6 n_step])';
[jac_dro] = Jacobi_Calc(x0_dro(1:n_int:n_step,:),mu);

% Define a colormap
cmap = colormap(jet(n_str));

% Sort the Jacobi constant for colormap (so colors depend on the Jacobi constant of a particular curve)
[jac_sorted,jac_sorted_ind] = sort(jac_dro);

for ii = 1:n_int:n_step
    
    % Iterate iteration counter
    iter = iter + 1;
    
    % Extract design variable vector for current step
    X0_i = X_store{ii};
    T_half_i = Thalf_store{ii};
    
    % Integrate converged result guess
    [x_fin,t_init] = GslInteg_STM([X0_i phi0],[0 2*T_half_i],mu);
    
    % Plot converged multiple shooting result
%     plot(x_fin(1,1).*l_ch,x_fin(1,2).*l_ch,'.r');
%     plot(x_fin(end,1).*l_ch,x_fin(end,2).*l_ch,'.g');
%     plot(x_fin(:,1).*l_ch,x_fin(:,2).*l_ch,'b');
    cmap_ind = find(iter == jac_sorted_ind);
    plot(x_fin(:,1).*l_ch,x_fin(:,2).*l_ch,'color',cmap(cmap_ind,:));

    % Store orbit initial conditions and period
    x0_store(iter,:) = X0_i;
    tf_store(iter,:) = 2*T_half_i;
    
    % Assemble and store the STM
    STM_store(:,:,iter) = reshape(x_fin(end,7:42),[6 6])';
    
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
caxis([min(jac_dro) max(jac_dro)]);
h = colorbar;
ylabel(h,'Jacobi Constant [nd]')
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

% Inex to plot up to for orbit stability indices
stp = n_str;

% Plot Stability indices 
figure (2)

% Index 1
subplot(3,1,1)
hold on
grid on
% plot(r_peri_store(1:stp),stable_ind(1:stp,1),'.b','MarkerSize',10);
% plot(r_peri_store(1:stp),stable_ubnd(1:stp),':k','LineWidth',1);
% plot(r_peri_store(1:stp),stable_lbnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_ind(1:stp,1),'.b','MarkerSize',10);
plot(Jac,stable_ubnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_lbnd(1:stp),':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Perilune Radius [km]')
ylabel('|\nu_1|')
title('DRO - Stability Index, \nu, as a Function of JC')
hold off

% Index 2
subplot(3,1,2)
hold on
grid on
% plot(r_peri_store(1:stp),stable_ind(1:stp,2),'.b','MarkerSize',10);
% plot(r_peri_store(1:stp),stable_ubnd(1:stp),':k','LineWidth',1);
% plot(r_peri_store(1:stp),stable_lbnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_ind(1:stp,2),'.b','MarkerSize',10);
plot(Jac,stable_ubnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_lbnd(1:stp),':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Perilune Radius [km]')
ylabel('|\nu_2|')
hold off

% Index 3
subplot(3,1,3)
hold on
grid on
% plot(r_peri_store(1:stp),stable_ind(1:stp,3),'.b','MarkerSize',10);
% plot(r_peri_store(1:stp),stable_ubnd(1:stp),':k','LineWidth',1);
% plot(r_peri_store(1:stp),stable_lbnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_ind(1:stp,3),'.b','MarkerSize',10);
plot(Jac,stable_ubnd(1:stp),':k','LineWidth',1);
plot(Jac,stable_lbnd(1:stp),':k','LineWidth',1);
xlabel('Jacobi Constant [nd]')
% xlabel('Perilune Radius [km]')
ylabel('|\nu_3|')
hold off

% Plot Jacobi constant as a function of perilune radius
figure (3)
hold on
grid on
plot(r_peri_store,Jac,'.b','MarkerSize',10);
xlabel('Perilune Radius [km]')
ylabel('Jacobi Constant [nd]')
title('DRO - Jacobi Constant as a Function of Perilune Radius')
hold off

% Plot Jacobi constant as a function of orbtial period
figure (4)
hold on
grid on
plot(tf_store.*t_ch./86400,Jac,'.b','MarkerSize',10);
xlabel('Orbital Period [days]')
ylabel('Jacobi Constant [nd]')
title('DRO - Jacobi Constant as a Function of Orbital Period')
hold off

% Plot Orbital Period as a function of Jacobi Constant
figure (5)
hold on
grid on
plot(Jac,tf_store.*t_ch./86400,'.b','MarkerSize',10);
xlabel('Jacobi Constant [nd]')
ylabel('Orbital Period [days]')
title('DRO - Jacobi Constant as a Function of Orbital Period')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Periodic Orbit Family as .mat File  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('EM_DRO_500_orbits','x0_store','tf_store')