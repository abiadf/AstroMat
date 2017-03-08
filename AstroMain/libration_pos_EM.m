function [L] = libration_pos_EM(mu)
% function [L] = libration_pos(mu)
% 
% Given the CR3BP mass ratio mu this code calculates the 3D position 
% coordinates of all five libration points in the Earth-Moon system. 
%
% INPUTS:
%    mu   	CR3BP Mass Ratio
%
% OUTPUTS:
%    L      Matrix containing position coordinates of libration points (5x3)
%
% Written by R. Pritchett, 7/24/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 Position %%

%Compose gam1 first guess for each system
gam1_guess = .25;
%Newton-Raphson Iteration
iter_1 = 0;
delt_gam_1 = 1;
gam1_new = gam1_guess;
while delt_gam_1 > 10^-12
    gam1 = gam1_new;
    num_1 = ((mu-1)/(1-gam1)^2)+(mu/gam1^2)+1-mu-gam1;
    den_1 = (2*(mu-1)/(1-gam1)^3)-((2*mu)/(gam1^3))-1;
    gam1_new = gam1-(num_1/den_1);
    delt_gam_1 = abs(gam1_new-gam1);
    iter_1 = iter_1+1;
end
gam1_fin = gam1_new; %final nondimensional gam1 values
%Calculate nondimensional L1 position
L1 = 1-mu-gam1_fin; %final nondimensional L1 values

%% L2 Position %%

%Set gam2 first guess
gam2_guess = 0.25;
%Newton-Raphson Iteration
iter_2 = 0;
delt_gam_2 = 1;
gam2_new = gam2_guess;
while delt_gam_2 > 10^-12
    gam2 = gam2_new;
    num_2 = ((mu-1)/(1+gam2)^2)-(mu/gam2^2)+1-mu+gam2;
    den_2 = (2*(1-mu)/(1+gam2)^3)+((2*mu)/(gam2^3))+1;
    gam2_new = gam2-(num_2/den_2);
    delt_gam_2 = abs(gam2_new-gam2);
    iter_2 = iter_2+1;
end
gam2_fin = gam2_new; %final nondimensional gam2 values
%Calculate nondimensional L2 position
L2 = 1-mu+gam2_fin; %final nondimensional L2 values

%% L3 Position %%

%Compose gam3 first guess for each system
gam3_guess = .25;
%Newton-Raphson Iteration
iter_3 = 0;
delt_gam_3 = 1;
gam3_new = gam3_guess;
while delt_gam_3 > 10^-12
    gam3 = gam3_new;
    num_3 = ((1-mu)/gam3^2)+(mu/(1+gam3)^2)-mu-gam3;
    den_3 = (-2*(1-mu)/(gam3)^3)-((2*mu)/((1+gam3)^3))-1;
    gam3_new = gam3-(num_3/den_3);
    delt_gam_3 = abs(gam3_new-gam3);
    iter_3 = iter_3+1;
end
gam3_fin = gam3_new; %final nondimensional gam3 values
%Calculate nondimensional L3 position
L3 = -mu-gam3_fin; %final nondimensional L3 values

%% Calculate L4 Position

L4_x = .5-mu;
L4_y = sqrt(3)/2; % note this is only true for the Earth-Moon system

%% Calculate L5 Position

L5_x = .5-mu;
L5_y = -sqrt(3)/2; % note this is only true for the Earth-Moon system

%% Output
L = [L1 0 0;L2 0 0;L3 0 0;L4_x L4_y 0;L5_x L5_y 0];