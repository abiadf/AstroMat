function [X,T_half,DF_fin,sshoot] = SingShoot_Lyap_NatParam(X0,T0_half,sshoot)
% function [X,T_half,DF_fin,sshoot] = SingShoot_Lyap_NatParam(X0,T0_half,sshoot)
% 
% Given a vector of initial conditions, a half period time and a structure 
% of parameters this function solves a single shooting problem for natural 
% motion in the CR3BP. Specifically, this algorithm targets periodic orbits
% that are planar and symmetric about the x-axis, this includes all three 
% families of Lyapunov orbits as well as the DRO's.
%
% INPUTS:
%    X0          initial guess at initial conditions of desired orbit 
%    T0_half     initial guess at half period of desired orbit
%    sshoot      structure containing single shooting parameters
%
% OUTPUTS:
%    X0          initial conditions of desired orbit 
%    T0_half     half period of desired orbit
%    DF_fin      DF matrix used for final update needed for use in Pseudo-Arclength Continuation method
%    sshoot      structure containing single shooting parameters
%
% Written by R. Pritchett, 02/13/17
% Last Update: R. Pritchett, 02/13/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from sshoot stucture
newt_tol = sshoot.newt_tol;
atten_tol = sshoot.atten_tol;
atten = sshoot.atten;
n_state = sshoot.n_state;
phi0 = sshoot.phi0;
mu = sshoot.mu;

% Initialize Newton's method iteration counter
newti = 0;

% Print header for output table
fprintf('\n                                                                                   \n')
fprintf(' len(X)  len(F)    max(|dX|)          max(|F|)      Iter  \n')
 
% Define design variable vector
Z = [X0(5); T0_half];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Form Initial Constraint Vector and Calculate Partials %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup Integration
tspan = [0 Z(2)]; 
ic = [X0(1:4), Z(1), X0(6), phi0];
event_id = 2; % row vector that identifies the event type
event_state = 0; % row vector defining desired event states
dir = 0; % row vector indicating desired event crossing directions (+1 - positive crossing, -1 - negative crossing, 0 - all crossings).
stop = 1; % row vector indicating whether integration should terminate at the events

% Integrate Trajectory
[x,~,~,~,~] = GslInteg_STM_events(ic,tspan,mu,event_id,event_state,dir,stop);
xf = x(end,:); % states at final time
phif = x(end,n_state+1:end); % STM at final time

% Calculate constraint vector
F = [xf(2); xf(4)]; % y and xdot at final time must equal zero

% Calculate error
error = max(abs(F));

% Calculate length of design variable and constraint vectors
Zl = length(Z);
Fl = length(F);

% Form STM
phif = reshape(phif,[n_state, n_state])'; % STM at final time

% Calculate time derivatives at final time
[xfdot] = EOM_CR3BP(0,xf,mu); % Note: input value of t is arbitrary

% Form Jacobian from STM and time derivatives
DF_new = [phif(2,5), xfdot(2);...
          phif(4,5), xfdot(4)];
  
% Use Newton's Method to calculate update
dZ = -DF_new'*((DF_new*DF_new')\F);
maxdZ = max(abs(dZ)); % calculate max change in X for print output

newti = newti + 1; % iterate counter by one

%Print output for one iteration of Newton's method
fprintf(' %3.0f  %6.0f  %16.8e  %16.8e  %5.0f \n',Zl,Fl,maxdZ,error,newti)    

% Reset DF_new to DF (placed here in case Newton's Method while loop is not entered)
DF = DF_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Newton's Method %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while error > newt_tol
    
    % Reset DF_sparse_new to DF_sparse
    DF = DF_new;
    
    if max(dZ) > atten_tol
        Z = Z + (1/atten).*dZ; % take quarter step
    else
        Z = Z + dZ; % take a full step
    end

    % Setup Integration
    tspan = [0 Z(2)]; 
    ic = [X0(1:4), Z(1), X0(6), phi0];
    event_id = 2; % row vector that identifies the event type
    event_state = 0; % row vector defining desired event states
    dir = 0; % row vector indicating desired event crossing directions (+1 - positive crossing, -1 - negative crossing, 0 - all crossings).
    stop = 1; % row vector indicating whether integration should terminate at the events

    % Integrate Trajectory
    [x,~,~,~,~] = GslInteg_STM_events(ic,tspan,mu,event_id,event_state,dir,stop);
    xf = x(end,:); % states at final time
    phif = x(end,n_state+1:end); % STM at final time

    % Calculate constraint vector
    F = [xf(2); xf(4)]; % y and xdot at final time must equal zero

    % Calculate error
    error = max(abs(F));

    % Calculate length of design variable and constraint vectors
    Zl = length(Z);
    Fl = length(F);

    % Form STM
    phif = reshape(phif,[n_state, n_state])'; % STM at final time

    % Calculate time derivatives at final time
    [xfdot] = EOM_CR3BP(0,xf,mu); % Note: input value of t is arbitrary

    % Form Jacobian from STM and time derivatives
    DF_new = [phif(2,5), xfdot(2);...
              phif(4,5), xfdot(4)];

    % Use Newton's Method to calculate update
    dZ = -DF_new'*((DF_new*DF_new')\F);
    maxdZ = max(abs(dZ)); % calculate max change in X for print output

    newti = newti + 1; % iterate counter by one

    %Print output for one iteration of Newton's method
    fprintf(' %3.0f  %6.0f  %16.8e  %16.8e  %5.0f \n',Zl,Fl,maxdZ,error,newti)    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Output IC and Half Period %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define output orbit IC
X = [X0(1:4) Z(1) X0(6)];

% Define output orbit half period
T_half = Z(2);

% Save DF matrix used for final update for use in Pseudo-Arclength Continuation method
DF_fin = DF;

% Save number of Newton's method iterations required to converge
sshoot.newti = newti;
