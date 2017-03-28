% Plot_Example_Direct_DRO2NRHO_EM.m
%
% This script plots the results of the associated script 
% Example_Direct_DRO2NRHO_EM.m which employs collocation and direct 
% optimization to compute low-thrust transfers between DRO and L1 and L2
% NRHO orbits
%
% Originally Written by: R. Pritchett, 09/23/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Configuration Space Plot %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot Setup %%%
figure(1)
hold on
grid on
axis equal

% Define initial view
% view(18,16)
view([-32,18])
% view(0,90) % XY Plane View
% view(0,0) % XZ View
% view(90,0) % YZ Plane View

% Plot Earth - realistic image, function from Bonnie Prado
bodyplot('Textures\','Earth',r_E,-mu.*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Moon - realistic image, function from Bonnie Prado
bodyplot('Textures\','Moon',r_M,(1-mu).*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Libration Points
dx = 0.02; % text offset in x
dy = 0.02; % text offset in y
dz = 0.0; % text offset in z
plot3(L(1,1).*l_ch,L(1,2).*l_ch,L(1,3).*l_ch,'sk') % plot L1
text((L(1,1)-dx).*l_ch,(L(1,2)+dy).*l_ch,(L(1,3)+dz).*l_ch,'L_1') % plot L1
plot3(L(2,1).*l_ch,L(2,2).*l_ch,L(2,3).*l_ch,'sk') % plot L2
text((L(2,1)-dx).*l_ch,(L(2,2)+dy).*l_ch,(L(2,3)+dz).*l_ch,'L_2') % plot L2
% plot3(L(3,1).*l_ch,L(3,2).*l_ch,L(3,3).*l_ch,'sk') % plot L3
% text((L(3,1)+dx).*l_ch,(L(3,2)+dy).*l_ch,(L(3,3)+dz).*l_ch,'L_3') % plot L3
% plot3(L(4,1).*l_ch,L(4,2).*l_ch,L(4,3).*l_ch,'sk') % plot L4
% text((L(4,1)+dx).*l_ch,(L(4,2)+dy).*l_ch,(L(4,3)+dz).*l_ch,'L_4') % plot L4
% plot3(L(5,1).*l_ch,L(5,2).*l_ch,L(5,3).*l_ch,'sk') % plot L5
% text((L(5,1)+dx).*l_ch,(L(5,2)+dy).*l_ch,(L(5,3)+dz).*l_ch,'L_5') % plot L5

% Add labels to plot
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')

% Plot Periodic Orbits
plot3(Xprop_init(:,1).*l_ch,Xprop_init(:,2).*l_ch,Xprop_init(:,3).*l_ch,':k','LineWidth',1); % initial halo orbit
plot3(Xprop_fin(:,1).*l_ch,Xprop_fin(:,2).*l_ch,Xprop_fin(:,3).*l_ch,':k','LineWidth',1); % final halo orbit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Converged Solution in Configuration Space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process direct transcription results
[ZPlot_fin] = Z2Plot(Z,x_bnd,t_bnd,t_var,colt);

% Extract data from Z2Plot structure for plotting
x_bnd_plot_fin = ZPlot_fin.x_bnd_plot;
x_var_plot_fin = ZPlot_fin.x_var_plot;
u_var_fin = ZPlot_fin.u_var;
T_var_fin = ZPlot_fin.T_var;
uT_var_fin = ZPlot_fin.uT_var;
pos_bnd_dim_plot = ZPlot_fin.pos_bnd_dim_plot;
vel_bnd_dim_plot = ZPlot_fin.vel_bnd_dim_plot;
mass_bnd_dim_plot = ZPlot_fin.mass_bnd_dim_plot;
pos_var_dim_plot = ZPlot_fin.pos_var_dim_plot;
vel_var_dim_plot = ZPlot_fin.vel_var_dim_plot;
mass_var_dim_plot = ZPlot_fin.mass_var_dim_plot;
T_var_fin_dim = ZPlot_fin.T_var_dim; 
uT_var_fin_dim = ZPlot_fin.uT_var_dim; 
t_var_hist = ZPlot_fin.t_var_hist;
T_var_hist = ZPlot_fin.T_var_hist;
uT_var_hist = ZPlot_fin.uT_var_hist;

% % Plot boundary nodes
% plot3(x_bnd_plot_fin(:,1).*l_ch,x_bnd_plot_fin(:,2).*l_ch,x_bnd_plot_fin(:,3).*l_ch,'+k','MarkerSize',8);

% Plot initial and final boundary nodes only
plot3(x_bnd_plot_fin(1,1).*l_ch,x_bnd_plot_fin(1,2).*l_ch,x_bnd_plot_fin(1,3).*l_ch,'.k','MarkerSize',12);
plot3(x_bnd_plot_fin(end,1).*l_ch,x_bnd_plot_fin(end,2).*l_ch,x_bnd_plot_fin(end,3).*l_ch,'.k','MarkerSize',12);

% Plot variable nodes
% plot3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,'.b','MarkerSize',8);
% plot3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,'b','LineWidth',2);

% % Plot thrust vectors
% quiver3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,...
%     uT_var_fin(:,1),uT_var_fin(:,2),uT_var_fin(:,3),'r'); % from indirect optimization result

% Plot results as thrust or coast arcs
x_traj = ZPlot_fin.x_traj;
t_traj = ZPlot_fin.t_traj;
traj_TorC = ZPlot_fin.traj_TorC;

% Define logical variable indicating whether coast or thrust period occurs
nocst = true;
nothrst = true;

for ii = 1:colt.n_seg
    
    % Extract data for current trajectory segment
    x_traj_i = x_traj{ii};
    traj_TorC_i = traj_TorC{ii};
    
    % Determine whether thrust or coast arc and plot
    scl = 100000;
    if strcmp(traj_TorC_i,'T')
        h1 = plot3(x_traj_i(:,1).*l_ch,x_traj_i(:,2).*l_ch,x_traj_i(:,3).*l_ch,'r','LineWidth',1);
        ind_arrow = floor(size(x_traj_i,1)/2);
        quiver3(x_traj_i(ind_arrow,1).*l_ch,x_traj_i(ind_arrow,2).*l_ch,x_traj_i(ind_arrow,3).*l_ch,...
            x_traj_i(ind_arrow,4).*scl,x_traj_i(ind_arrow,5).*scl,x_traj_i(ind_arrow,6).*scl,'r')
        nothrst = false;
    else
        h2 = plot3(x_traj_i(:,1).*l_ch,x_traj_i(:,2).*l_ch,x_traj_i(:,3).*l_ch,'b','LineWidth',1);
%         ind_arrow = floor(size(x_traj_i,1)/2);
%         quiver3(x_traj_i(ind_arrow,1).*l_ch,x_traj_i(ind_arrow,2).*l_ch,x_traj_i(ind_arrow,3).*l_ch,...
%             x_traj_i(ind_arrow,4).*scl,x_traj_i(ind_arrow,5).*scl,x_traj_i(ind_arrow,6).*scl,'b') 
        nocst = false;
    end
end

% Add legend
if nocst
    legend(h1,{'Thrust'})
elseif nothrst
    legend(h2,{'Coast'})
else
    legend([h1 h2],{'Thrust','Coast'})
end

hold off
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Thrust and Specific Impulse %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (2)

% Thrust magnitude
% subplot(1,2,1)
hold on
grid on
h1 = plot(t_var_hist.*(t_ch/86400),T_var_hist,'k');
h2 = plot(t_var_hist.*(t_ch/86400),uT_var_hist(1,:),'b');
h3 = plot(t_var_hist.*(t_ch/86400),uT_var_hist(2,:),'r');
h4 = plot(t_var_hist.*(t_ch/86400),uT_var_hist(3,:),'g');
legend([h1,h2,h3,h4],{'|T|','T1','T2','T3'})
xlabel('Time [days]')
ylabel('Thrust [N]')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Approximate "Delta-V" Associated with Transfer %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltV_approx = colt.ce_dim*log(mass_bnd_dim_plot(1)/mass_bnd_dim_plot(end)); % [m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Package Data for Export to Other Users %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rename variables for export
X_ppt = x_bnd_plot_fin;
T_ppt = t_bnd;

% save('DRO2L1NRHO_OptwColl_rev11_Tmax10_n60_Export_BndPts','mu','X_init',...
%     'T_init','X_fin','T_fin','X_ppt','T_ppt')
