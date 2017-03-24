% Plot_Example_Direct_H2HMan_EM.m
%
% This script plots the results of the associated script 
% Example_Direct_H2HMan_EM.m which employs collocation and direct 
% optimization to compute low-thrust transfers between periodic halo orbits
% that leverage invariant manifold structures.
%
% Originally Written by: R. Pritchett, 09/23/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data Necessary for Plotting from Structures %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Matrices from integration timestep structure
% XT = intstp_dscrt.XT;
% lambda_intstp = intstp_dscrt.lambda_intstp;
% lambda_m_intstp = intstp_dscrt.lambda_m_intstp;
% t_prop = intstp_dscrt.t_prop;
% T_intstp = intstp_dscrt.T_intstp; 
% T_intstp_dim = intstp_dscrt.T_intstp_dim; 
% u_intstp = intstp_dscrt.u_intstp; 
% P_intstp = intstp_dscrt.P_intstp;
% 
% % Matrices from boundary node structure
% XT_bnd = bndnd_dscrt.XT_bnd;
% lambda_bnd = bndnd_dscrt.lambda_bnd;
% lambda_m_bnd = bndnd_dscrt.lambda_m_bnd; 
% t_prop_bnd = bndnd_dscrt.t_bnd;
% T_bnd = bndnd_dscrt.T_bnd; 
% T_bnd_dim = bndnd_dscrt.T_bnd_dim; 
% u_bnd = bndnd_dscrt.u_bnd; 
% P_bnd = bndnd_dscrt.P_bnd;
% 
% % Matrices from variable node structure
% XT_var = varnd_dscrt.XT_var;
% lambda_var = varnd_dscrt.lambda_var;
% lambda_m_var = varnd_dscrt.lambda_m_var; 
% t_prop_var = varnd_dscrt.t_var;
% T_var = varnd_dscrt.T_var; 
% T_var_dim = varnd_dscrt.T_var_dim; 
% u_var = varnd_dscrt.u_var; 
% P_var = varnd_dscrt.P_var;

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
view(0,90)

% Plot Moon - realistic image, function from Bonnie Prado
bodyplot('Textures\','Moon',r_M,(1-mu).*l_ch,0,0,0.9,[1 0 0],0); % create 3D surface from sphere coordinates

% Plot Libration Points
dx = +0.03; % text offset in x
dy = 0.02; % text offset in y
dz = 0.0; % text offset in z
plot3(L(1,1).*l_ch,L(1,2).*l_ch,L(1,3).*l_ch,'sk') % plot L1
text((L(1,1)+0.3*dx).*l_ch,(L(1,2)+dy).*l_ch,(L(1,3)+dz).*l_ch,'L_1') % plot L1
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
%% Plot Discretized Initial Guess in Configuration Space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Process direct transcription results
% [ZPlot_init] = Z2Plot(Z0,x_bnd,t_bnd,t_var,colt);
% 
% % Save results in structure for output to plotting script
% x_bnd_plot_init = ZPlot_init.x_bnd_plot;
% x_var_plot_init = ZPlot_init.x_var_plot;
% u_var_init = ZPlot_init.u_var;
% T_var_init = ZPlot_init.T_var;
% P_var_init = ZPlot_init.P_var;
% uT_var_init = ZPlot_init.uT_var;
% X_tau1_init = ZPlot_init.X_tau1;
% X_tau2_init = ZPlot_init.X_tau2;
% X_alph1_init = ZPlot_init.X_alph1;
% X_alph2_init = ZPlot_init.X_alph2;

% % Plot coast arcs
% plot3(X_tau1_init(:,1),X_tau1_init(:,2),X_tau1_init(:,3),':b','LineWidth',2) % orbit 1 coast arc
% plot3(X_tau1_init(end,1),X_tau1_init(end,2),X_tau1_init(end,3),'.b','MarkerSize',20) % end of orbit 1 coast arc
% plot3(X_tau2_init(:,1),X_tau2_init(:,2),X_tau2_init(:,3),':r','LineWidth',2) % orbit 2 coast arc
% plot3(X_tau2_init(end,1),X_tau2_init(end,2),X_tau2_init(end,3),'.r','MarkerSize',20) % end of orbit 2 coast arc

% % Plot manifold arcs
% plot3(X_alph1_init(:,1),X_alph1_init(:,2),X_alph1_init(:,3),':c','LineWidth',2) % orbit 1 coast arc
% plot3(X_alph1_init(end,1),X_alph1_init(end,2),X_alph1_init(end,3),'.c','MarkerSize',20) % end of orbit 1 coast arc
% plot3(X_alph2_init(:,1),X_alph2_init(:,2),X_alph2_init(:,3),':m','LineWidth',2) % orbit 2 coast arc
% plot3(X_alph2_init(end,1),X_alph2_init(end,2),X_alph2_init(end,3),'.m','MarkerSize',20) % end of orbit 2 coast arc

% Plot boundary nodes
% plot3(x_bnd_plot_init(:,1),x_bnd_plot_init(:,2),x_bnd_plot_init(:,3),'ok','MarkerSize',5);

% Plot variable nodes
% plot3(x_var_plot_init(:,1),x_var_plot_init(:,2),x_var_plot_init(:,3),'or','MarkerSize',5);

% Plot thrust vectors
% quiver3(x_var_plot_init(:,1),x_var_plot_init(:,2),x_var_plot_init(:,3),...
%     uT_var_init(:,1),uT_var_init(:,2),uT_var_init(:,3),'g'); % from indirect optimization result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Converged Solution in Configuration Space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process direct transcription results
[ZPlot_fin] = Z2Plot(Z,x_bnd,t_bnd,t_var,colt);

% Save results in structure for output to plotting script
x_bnd_plot_fin = ZPlot_fin.x_bnd_plot;
x_var_plot_fin = ZPlot_fin.x_var_plot;
u_var_fin = ZPlot_fin.u_var;
T_var_fin = ZPlot_fin.T_var;
uT_var_fin = ZPlot_fin.uT_var;
% X_tau1_fin = ZPlot_fin.X_tau1;
% X_tau2_fin = ZPlot_fin.X_tau2;
% X_alph1_fin = ZPlot_fin.X_alph1;
% X_alph2_fin = ZPlot_fin.X_alph2;
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

% % Plot coast arcs
% plot3(X_tau1_fin(:,1),X_tau1_fin(:,2),X_tau1_fin(:,3),':b','LineWidth',2) % orbit 1 coast arc
% plot3(X_tau1_fin(end,1),X_tau1_fin(end,2),X_tau1_fin(end,3),'.b','MarkerSize',20) % end of orbit 1 coast arc
% plot3(X_tau2_fin(:,1),X_tau2_fin(:,2),X_tau2_fin(:,3),':r','LineWidth',2) % orbit 2 coast arc
% plot3(X_tau2_fin(end,1),X_tau2_fin(end,2),X_tau2_fin(end,3),'.r','MarkerSize',20) % end of orbit 2 coast arc

% % Plot manifold arcs
% plot3(X_alph1_fin(:,1),X_alph1_fin(:,2),X_alph1_fin(:,3),':c','LineWidth',2) % orbit 1 coast arc
% plot3(X_alph1_fin(end,1),X_alph1_fin(end,2),X_alph1_fin(end,3),'.c','MarkerSize',20) % end of orbit 1 coast arc
% plot3(X_alph2_fin(:,1),X_alph2_fin(:,2),X_alph2_fin(:,3),':m','LineWidth',2) % orbit 2 coast arc
% plot3(X_alph2_fin(end,1),X_alph2_fin(end,2),X_alph2_fin(end,3),'.m','MarkerSize',20) % end of orbit 2 coast arc

% % Plot boundary nodes
% plot3(x_bnd_plot_fin(:,1).*l_ch,x_bnd_plot_fin(:,2).*l_ch,x_bnd_plot_fin(:,3).*l_ch,'+k','MarkerSize',8);

% Plot initial and final boundary nodes only
plot3(x_bnd_plot_fin(1,1).*l_ch,x_bnd_plot_fin(1,2).*l_ch,x_bnd_plot_fin(1,3).*l_ch,'.k','MarkerSize',12);
plot3(x_bnd_plot_fin(end,1).*l_ch,x_bnd_plot_fin(end,2).*l_ch,x_bnd_plot_fin(end,3).*l_ch,'.k','MarkerSize',12);

% Plot variable nodes
% plot3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,'.b','MarkerSize',8);
% plot3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,'b','LineWidth',2);

% Plot thrust vectors
quiver3(x_var_plot_fin(:,1).*l_ch,x_var_plot_fin(:,2).*l_ch,x_var_plot_fin(:,3).*l_ch,...
    uT_var_fin(:,1),uT_var_fin(:,2),uT_var_fin(:,3),'r'); % from indirect optimization result

% Plot results as thrust or coast arcs
x_traj = ZPlot_fin.x_traj;
t_traj = ZPlot_fin.t_traj;
traj_TorC = ZPlot_fin.traj_TorC;

for ii = 1:colt.n_seg
    
    % Extract data for current trajectory segment
    x_traj_i = x_traj{ii};
    traj_TorC_i = traj_TorC{ii};
    
    % Determine whether thrust or coast arc and plot
    if strcmp(traj_TorC_i,'T')
        h1 = plot3(x_traj_i(:,1).*l_ch,x_traj_i(:,2).*l_ch,x_traj_i(:,3).*l_ch,'r','LineWidth',1);
    else
        h2 = plot3(x_traj_i(:,1).*l_ch,x_traj_i(:,2).*l_ch,x_traj_i(:,3).*l_ch,'b','LineWidth',1);
    end
    
end

% Add legend, modify if no coast segments exist
if exist('h2','var')
    legend([h1 h2],{'Thrust','Coast'})
else
    legend(h1,{'Thrust'})
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




