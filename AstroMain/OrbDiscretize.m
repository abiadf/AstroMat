function [Z0,x_plot] = OrbDiscretize(X,T,mult)
% function [Z0,x_plot] = OrbDiscretize(X,T,mult)
%
% This function discretizes a trajectory into patch points and assembles a 
% vector of design variables, including the total trajectory time, that can 
% be passed as an initial guess to a multiple shooting method. 
%
% INPUTS:
%    X        vector of initial state values,(6x1)
%    T        orbital period
%    mult     structure containing multiple shooting parameters parameters      
%
% OUTPUTS:
%    Z0       column vector containing initial guess for a multiple 
%             shooting method.
%    x_plot   states at patch points collected in matrix form for simple plotting 
%
% Originally Written by: R. Pritchett, 02/09/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from mult stucture
n = mult.n;
n_state = mult.n_state;
options = mult.options;
mu = mult.mu;

% Define patch point times
t_ppt = linspace(0,T,n);

% Integrate orbit and record states at variable node times
[~,x_ppt] = CR3BP_EOM_NO_STM(X,t_ppt,mu);

% Reshape states into column vector of design variables
Z0 = reshape(x_ppt',[n*n_state 1]);

% Add total trajectory time to design variable vector
Z0 = [Z0;T];

% Reshape (if necessary) x_ppt into form amenable to plotting
x_plot = x_ppt;

