function [Xorb_traj,Torb,V_mono] = PeriodicProp(OrbIC,colt)
% function [Xorb_traj,Torb,V_mono] = PeriodicProp(OrbIC,phi0,options,mu)
%
% Given a set of initial position and velocity states as well as a multiple
% of the orbital period this function propagates a natural periodic orbit 
% and computes its monodromy matrix.
%
% INPUTS:
%    OrbIC      initial position, velocity, and period of a periodic orbit [7x1]   
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    Xorb_traj  matrix of position and velocity states at each time step along the periodic orbit
%    Torb       vector of time steps along the periodic orbit
%    V_mono     sorted monodromy matrix of the periodic orbit
%
% Originally Written by: R. Pritchett, 09/01/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
phi0 = colt.phi0; % initial value of the STM
options = colt.options; % numerical integration options
mu = colt.mu; % CR3BP mass ratio 

% Extract periodic orbit IC
Xorb0 = OrbIC(1:6);
Torb0 = OrbIC(7);

% Propagate orbit
[Xorb,Torb] = GslInteg_STM([Xorb0 phi0],[0 Torb0],mu);
Xorb_traj = Xorb(:,1:6); % extract position and velocity states of initial orbit
Xorb_STM = Xorb(:,7:42); % extract STM components of initial orbit

% Define monodromy matrix for orbit
mono_orb = reshape(Xorb_STM(end,:),[6 6])';

% Calculate eigenvalues and eigenvectors for monodromy matrices
[V_mono,~] = eig(mono_orb);
