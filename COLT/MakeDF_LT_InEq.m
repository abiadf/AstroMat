function [DF_InEq] = MakeDF_LT_InEq(colt)
% function [DF] = MakeDF_LT_InEq(Z,t_seg,t_seg_d,collmat,colt)
%
% This function calculates the partials of the inequality constraints for
% direct transcription problems involving low-thrust. Matrices of partials 
% in collocation problems are typically highly sparse, thus this function 
% calculates only the partials that are nonzero.
%
% INPUTS:
%    colt       structure containing collocation and optimization parameters
%
% OUTPUTS:
%    DF_InEq    nonzeros of the Jacobian matrix in column vector form
%
% Written by R. Pritchett, 02/02/2017
% Last Update: R. Pritchett, 02/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract necessary parameters from colt stucture
n_seg = colt.n_seg;

%Preallocate matrix of partial derivatives for inequality constraints
D_ineq = zeros(2,n_seg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Partials wrt Thrust Magnitude %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_ineq(1,:) = -ones(1,n_seg);
D_ineq(2,:) = ones(1,n_seg);
                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reshape components of DF into a single column vector 
DF_InEq =  reshape(D_ineq,[2*n_seg 1]);
      