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
mult.newt_tol = 1e-8; % tolerance for Newton's Method
mult.DiffType = 'Central'; % finite difference method for numerical jacobian
mult.n = 4; % number of patch points to use
mult.n_state = 6; % number of state variables
mult.n_slack = 0; % number of slack variables

%%% Boundary Constraints %%%

% Specify Boundary Constraint Case
mult.BndCase = 'FixEndPtOnly';
% mult.BndCase = 'Periodicity';

% Periodicity constraints
mult.period = [1 2 3 4 6]; % define which five states to include in periodicity constraint
mult.n_period = length(mult.period); % calculate number of periodicity constraints

% Initial Hyperplane Constraint
mult.i_hypr0 = [2 4 6]; % defines which component of the initial state is constrained to a hyperplane 
mult.x_hypr0 = [0 0 0]; % defines the value(s) of the component(s) that define the hyperplane 
mult.n_hypr0 = length(mult.i_hypr0); % calculate number of hyperplane constraints

% Final Hyperplane Constraint
mult.i_hyprf = [2 4 6]; % defines which component of the initial state is constrained to a hyperplane 
mult.x_hyprf = [0 0 0]; % defines the value(s) of the component(s) that define the hyperplane 
mult.n_hyprf = length(mult.i_hyprf); % calculate number of hyperplane constraints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define/Load Initial Guess and Discretize %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set of initial conditions from Tom Pavlak
load EM_L2_Lyap_250_orbits

% Define initial conditions
orb_num = 50;
X0 = u_0_save(orb_num,:);
T0 = t_f_save(orb_num);
tspan = [0 T0];

% Integrate initial guess
[t_init,x_init] = CR3BP_EOM_NO_STM(X0,tspan,mu);

% Discretize initial guess and assemble design variable vector
[Z0,x_plot] = OrbDiscretize(X0,T0,mult);

% Calculate length of design variable vector
mult.Zl = length(Z0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Guess %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
grid on
axis equal

plot(X0(1),X0(2),'*k');
plot(x_plot(:,1),x_plot(:,2),'.r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Converge Initial Solution and Calculate Null Space Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run multiple shooting algorithm
[Z,DF_fin] = MultShoot(Z0,mult);

% Save or load collocation results
% save('MSLyap_n10','Z','x_bnd','t_var','t_bnd','C')
% load MSLyap_n10.mat

% Convert DF_fin from sparse to full matrix
DF_fin_full = full(DF_fin);

% Calculate null vector from collocation results
delt_Z_prev = null(DF_fin_full);

% Define pseudo-arclength continuation scaling parameter
delt_s = .05;

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
n_step = 20;

% Preallocate cell arrays for storing results
Z_store = cell(1,n_step);

for ii = 1:n_step
 
    % Run multiple shooting algorithm
    [Z,DF_fin] = MultShoot_PArc(Z0,mult,parc);

    % Save or load collocation results
    % save('MSLyap_n10','Z','x_bnd','t_var','t_bnd','C')
    % load MSLyap_n10.mat

    % Store Results
    Z_store{ii} = Z;
    
    % Convert DF_fin from sparse to full matrix
    DF_fin_full = full(DF_fin);

    % Calculate null vector from collocation results
    delt_Z_prev = null(DF_fin_full);
    
    % Label Z from converged solution as Zprev, since new Z will be Zpls
    Z_prev = Z;
    
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

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;

% Convert column vector of design variables into matrix of state variables
x_ppt = reshape(Z(1:end-1),[n_state n]);

% Calculate patch point times using total time
T = Z(end);
t_ppt = linspace(0,T,n);

% Plot converged multiple shooting result
plot(x_ppt(1,:),x_ppt(2,:),'.b');
hold off
