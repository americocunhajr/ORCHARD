
% -----------------------------------------------------------------
%  post_orchard_ivp_kl.m
%
%  This script performs the post processing of simulation data
%  that comes from orchard nonlinear dynamics.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------


% plot left tire displacement
% ...........................................................
%gtitle = ' orchard sprayer track (left tire)';
gtitle = ' ';
xlab   = ' distance from origin (m)';
ylab   = ' displacement (m)';
xmin   = 0.0;
xmax   = vtrans*t1;
%ymin   = 0.0;
%ymax   = 0.2;
%xmin   = 'auto';
%xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_ye1_disp'];
flag   = 'eps';
fig0a  = graph_type1(tspan*vtrans,Ye1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig0a);
% ...........................................................


% plot right tire displacement
% ...........................................................
%gtitle = ' orchard sprayer track (right tire)';
gtitle = ' ';
xlab   = ' distance from origin (m)';
ylab   = ' displacement (m)';
xmin   = 0.0;
xmax   = vtrans*t1;
%ymin   = 0.0;
%ymax   = 0.2;
%xmin   = 'auto';
%xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_ye2_disp'];
flag   = 'eps';
fig0b  = graph_type1(tspan*vtrans,Ye2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig0b);
% ...........................................................


% plot trailer vertical displacement time series
% ...........................................................
%gtitle = ' trailer vertical displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = -0.1;
ymax   =  0.3;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_y1_disp'];
flag   = 'eps';
fig1a  = graph_type1(time,Qy1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................



% plot tower vertical displacement time series
% ...........................................................
%gtitle = ' tower vertical displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 2.5;
ymax   = 2.9;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_y2_disp'];
flag   = 'eps';
fig1e  = graph_type1(time,Qy2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1e);
% ...........................................................



% plot tower vertical displacement time series
% ...........................................................
%gtitle = ' tower vertical displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 2.5;
ymax   = 2.9;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_y12_disp'];
flag   = 'eps';
fig1e  = graph_type1(time,Qy2-Qy1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1e);
% ...........................................................


% plot trailer angular displacement time series
% ...........................................................
%gtitle = ' trailer angular displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = -0.2;
ymax   =  0.2;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_phi1_disp'];
flag   = 'eps';
fig1b  = graph_type1(time,Qphi1,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1b);
% ...........................................................




% plot tower angular displacement time series
% ...........................................................
%gtitle = ' tower angular displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = -0.2;
ymax   =  0.2;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_phi2_disp'];
flag   = 'eps';
fig1c  = graph_type1(time,Qphi2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1c);
% ...........................................................



% plot trailer angular displacement difference time series
% ...........................................................
%gtitle = ' trailer angular displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = -0.2;
ymax   =  0.2;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_phi12_disp'];
flag   = 'eps';
fig1cc  = graph_type1(time,Qphi2-Qphi1,gtitle,....
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1cc);
% ...........................................................




% plot tower horizontal displacement time series
% ...........................................................
%gtitle = ' tower horizontal displacement';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = -0.4;
ymax   =  0.4;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_x2_disp'];
flag   = 'eps';
fig1d  = graph_type1(time,Qx2,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1d);
% ...........................................................


% plot trailer vertical phase space trajectory
% ...........................................................
%gtitle = ' trailer vertical phase space trajectory';
gtitle = ' ';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = -0.1;
xmax   =  0.3;
ymin   = -2.0;
ymax   =  2.0;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_y1'];
flag   = 'eps';
fig2a  = graph_type1(Qy1(Nss:end),Qy1dot(Nss:end),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2a);
% ...........................................................


% plot trailer angular phase space trajectory
% ...........................................................
%gtitle = ' trailer angular phase space trajectory';
gtitle = ' ';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = -0.2;
xmax   =  0.2;
ymin   = -0.8;
ymax   =  0.8;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi1'];
flag   = 'eps';
fig2b  = graph_type1(Qphi1(Nss:end),Qphi1dot(Nss:end),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2b);
% ...........................................................


% plot tower angular phase space trajectory
% ...........................................................
%gtitle = ' tower angular phase space trajectory';
gtitle = ' ';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = -0.2;
xmax   =  0.2;
ymin   = -0.8;
ymax   =  0.8;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi2'];
flag   = 'eps';
fig2c  = graph_type1(Qphi2(Nss:end),Qphi2dot(Nss:end),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2c);
% ...........................................................


% plot tower horizontal phase space trajectory
% ...........................................................
%gtitle = ' tower horizontal phase space trajectory';
gtitle = ' ';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = -0.4;
xmax   =  0.4;
ymin   = -1.0;
ymax   =  1.0;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_x2'];
flag   = 'eps';
fig2d  = graph_type1(Qx2(Nss:end),Qx2dot(Nss:end),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2d);
% ...........................................................


% plot tower vertical phase space trajectory
% ...........................................................
%gtitle = ' tower vertical phase space trajectory';
gtitle = ' ';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = 2.5;
xmax   = 2.9;
ymin   = -2.0;
ymax   =  2.0;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_y2'];
flag   = 'eps';
fig2e  = graph_type1(Qy2(Nss:end),Qy2dot(Nss:end),gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2e);
% ...........................................................



% number of time steps to jump
Njump = 1;

% plot trailer vertical atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
%gtitle = ' trailer vertical atractor';
gtitle = ' ';
xmin   = -0.1;
xmax   =  0.3;
ymin   = -2.0;
ymax   =  2.0;
% xmin   = min(Qy1);
% xmax   = max(Qy1);
% ymin   = min(Qy1dot);
% ymax   = max(Qy1dot);
zmin   = min(time(Nss:end));
zmax   = max(time(Nss:end));
gname  = [num2str(case_name),'__atractor_y1'];
flag   = 'eps';
fig4a  = graph_3d_atractor(   Qy1(Nss:Njump:end),...
                          Qy1dot(Nss:Njump:end),...
                          time(Nss:Njump:end),...
                          time(Nss:Njump:end),...
                          gtitle,xlab,ylab,zlab,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4a);
% ...........................................................


% plot trailer angular atractor
% ...........................................................
xlab   = ' rotation (rad)';
ylab   = ' ang. velocity (rad/s)';
zlab   = ' time (s)';
%gtitle = ' trailer angular atractor';
gtitle = ' ';
xmin   = -0.2;
xmax   =  0.2;
ymin   = -0.8;
ymax   =  0.8;
% xmin   = min(Qphi1);
% xmax   = max(Qphi1);
% ymin   = min(Qphi1dot);
% ymax   = max(Qphi1dot);
zmin   = min(time(Nss:end));
zmax   = max(time(Nss:end));
gname  = [num2str(case_name),'__atractor_phi1'];
flag   = 'eps';
fig4b  = graph_3d_atractor(   Qphi1(Nss:Njump:end),...
                          Qphi1dot(Nss:Njump:end),...
                          time(Nss:Njump:end),...
                          time(Nss:Njump:end),...
                          gtitle,xlab,ylab,zlab,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4b);
% ...........................................................


% plot tower angular atractor
% ...........................................................
xlab   = ' rotation (rad)';
ylab   = ' ang. velocity (rad/s)';
zlab   = ' time (s)';
%gtitle = ' tower angular atractor';
gtitle = ' ';
xmin   = -0.2;
xmax   =  0.2;
ymin   = -0.8;
ymax   =  0.8;
% xmin   = min(Qphi2);
% xmax   = max(Qphi2);
% ymin   = min(Qphi2dot);
% ymax   = max(Qphi2dot);
zmin   = min(time(Nss:end));
zmax   = max(time(Nss:end));
gname  = [num2str(case_name),'__atractor_phi2'];
flag   = 'eps';
fig4c  = graph_3d_atractor(   Qphi2(Nss:Njump:end),...
                           Qphi2dot(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4c);
% ...........................................................


% plot tower horizontal atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
%gtitle = ' tower horizontal atractor';
gtitle = ' ';
xmin   = -0.4;
xmax   =  0.4;
ymin   = -1.0;
ymax   =  1.0;
% xmin   = min(Qx2);
% xmax   = max(Qx2);
% ymin   = min(Qx2dot);
% ymax   = max(Qx2dot);
zmin   = min(time(Nss:end));
zmax   = max(time(Nss:end));
gname  = [num2str(case_name),'__atractor_x2'];
flag   = 'eps';
fig4d  = graph_3d_atractor(   Qx2(Nss:Njump:end),...
                           Qx2dot(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4d);
% ...........................................................


% plot tower vertical atractor
% ...........................................................
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
zlab   = ' time (s)';
%gtitle = ' tower vertical atractor';
gtitle = ' ';
xmin   = 2.5;
xmax   = 2.9;
ymin   = -2.0;
ymax   =  2.0;
% xmin   = min(Qy2);
% xmax   = max(Qy2);
% ymin   = min(Qy2dot);
% ymax   = max(Qy2dot);
zmin   = min(time(Nss:end));
zmax   = max(time(Nss:end));
gname  = [num2str(case_name),'__atractor_y2'];
flag   = 'eps';
fig4e  = graph_3d_atractor(   Qy2(Nss:Njump:end),...
                           Qy2dot(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           time(Nss:Njump:end),...
                           gtitle,xlab,ylab,zlab,...
                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4e);
% ...........................................................




%time_8kmph = time;
%Qx2_8kmph  = Qx2;
%save('x2_corr_1m_v_08kmph.mat','time_8kmph','Qx2_8kmph')
%time_12kmph = time;
%Qx2_12kmph  = Qx2;
%save('x2_corr_1m_v_12kmph.mat','time_12kmph','Qx2_12kmph')
%time_16kmph = time;
%Qx2_16kmph  = Qx2;
%save('x2_corr_1m_v_16kmph.mat','time_16kmph','Qx2_16kmph')


% plot tower horizontal displacement time series
% ...........................................................
%gtitle = ' tower horizontal displacement';
gtitle = ' ';
leg1 = '08km/h';
leg2 = '12km/h';
leg3 = '16km/h';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = -0.4;
ymax   =  0.4;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_x2_disp_velocity'];
flag   = 'eps';
fig5a  = graph_type3( time_8kmph, Qx2_8kmph,...
                     time_12kmph,Qx2_12kmph,...
                     time_16kmph,Qx2_16kmph,...
                     gtitle,leg1,leg2,leg3,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5a);
% ...........................................................




time_01m = time;
Qx2_01m  = Qx2;
save('x2_corr_01m_v_12kmph.mat','time_01m','Qx2_01m')
time_1m = time;
Qx2_1m  = Qx2;
save('x2_corr_1m_v_12kmph.mat','time_1m','Qx2_1m')
time_10m = time;
Qx2_10m  = Qx2;
save('x2_corr_10m_v_12kmph.mat','time_10m','Qx2_10m')
time_50m = time;
Qx2_50m  = Qx2;
save('x2_corr_50m_v_12kmph.mat','time_50m','Qx2_50m')


% plot tower horizontal displacement time series
% ...........................................................
%gtitle = ' tower horizontal displacement';
gtitle = ' ';
leg1 = '0.1m';
leg2 = '  1m';
leg3 = ' 10m';
leg4 = ' 50m';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = -0.4;
ymax   =  0.4;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__tseries_x2_disp_correlation'];
flag   = 'eps';
fig5b  = graph_type4(time_01m,Qx2_01m,...
                      time_1m,Qx2_1m,...
                     time_10m,Qx2_10m,...
                     time_50m,Qx2_50m,...
                     gtitle,leg1,leg2,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5b);
% ...........................................................

