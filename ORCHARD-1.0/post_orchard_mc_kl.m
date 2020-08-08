
% -----------------------------------------------------------------
%  post_orchard_mc_kl.m
%
%  This script performs the post processing of simulation data
%  that comes from orchard nonlinear stochastic dynamics.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 15, 2017
% -----------------------------------------------------------------


samp_pts = randi([1 Ns],1,5);


% plot trailer vertical displacement confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' trailer vertical displacement';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_y1_disp'];
gname1 = [num2str(case_name),'__st_y1_disp'];
flag   = 'eps';
% fig1a  = graph_ci1(time,y1_smp_avg,...
%                    y1_upp,...
%                    y1_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1a  = graph_ci1N(time,MC_y1(samp_pts,:),...
                    y1_upp,...
                    y1_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
                
fig1aa = graph_type2(time,y1_smp_avg,...
                     time,y1_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1a);
% -----------------------------------------------------------


% plot trailer angular displacement confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' trailer angular displacement';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_phi1_disp'];
gname1 = [num2str(case_name),'__st_phi1_disp'];
flag   = 'eps';
% fig1b  = graph_ci1(time,phi1_smp_avg,...
%                    phi1_upp,...
%                    phi1_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1b  = graph_ci1N(time,MC_phi1(samp_pts,:),...
                    phi1_upp,...
                    phi1_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1bb = graph_type2(time,phi1_smp_avg,...
                     time,phi1_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1b);
% -----------------------------------------------------------



% plot tower angular displacement confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower angular displacement';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' rotation (rad)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_phi2_disp'];
gname1 = [num2str(case_name),'__st_phi2_disp'];
flag   = 'eps';
% fig1c  = graph_ci1(time,phi2_smp_avg,...
%                    phi2_upp,...
%                    phi2_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1c  = graph_ci1N(time,MC_phi2(samp_pts,:),...
                     phi2_upp,...
                     phi2_low,...
                     gtitle,leg1,leg2,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1cc = graph_type2(time,phi2_smp_avg,...
                     time,phi2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1c);
% -----------------------------------------------------------


% plot tower horizontal displacement confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower horizontal displacement';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_x2_disp'];
gname1 = [num2str(case_name),'__st_x2_disp'];
flag   = 'eps';
% fig1d  = graph_ci1(time,x2_smp_avg,...
%                    x2_upp,...
%                    x2_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1d  = graph_ci1N(time,MC_x2(samp_pts,:),...
                    x2_upp,...
                    x2_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1dd = graph_type2(time,x2_smp_avg,...
                     time,x2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1d);
% -----------------------------------------------------------



% plot tower vertical displacement confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower vertical displacement';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_y2_disp'];
gname1 = [num2str(case_name),'__st_y2_disp'];
flag   = 'eps';
% fig1d  = graph_ci1(time,y2_smp_avg,...
%                    y2_upp,...
%                    y2_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1e  = graph_ci1N(time,MC_y2(samp_pts,:),...
                    y2_upp,...
                    y2_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1ee = graph_type2(time,y2_smp_avg,...
                     time,y2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1e);
% -----------------------------------------------------------


% plot trailer vertical velocity confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' trailer vertical velocity';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' velocity (m/s)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_y1_velo'];
gname1 = [num2str(case_name),'__st_y1_velo'];
flag   = 'eps';
% fig2a  = graph_ci1(time,y1dot_smp_avg,...
%                    y1dot_upp,...
%                    y1dot_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2a  = graph_ci1N(time,MC_y1dot(samp_pts,:),...
                    y1dot_upp,...
                    y1dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2aa = graph_type2(time,y1dot_smp_avg,...
                     time,y1dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2a);
% -----------------------------------------------------------



% plot trailer angular velocity confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' trailer angular velocity';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' angular velocity (rad/s)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_phi1_velo'];
gname1 = [num2str(case_name),'__st_phi1_velo'];
flag   = 'eps';
% fig2b  = graph_ci1(time,phi1dot_smp_avg,...
%                    phi1dot_upp,...
%                    phi1dot_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2b  = graph_ci1N(time,MC_phi1dot(samp_pts,:),...
                    phi1dot_upp,...
                    phi1dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2bb = graph_type2(time,phi1dot_smp_avg,...
                     time,phi1dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2b);
% -----------------------------------------------------------



% plot tower angular velocity confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower angular velocity';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' angular velocity (rad/s)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_phi2_velo'];
gname1 = [num2str(case_name),'__st_phi2_velo'];
flag   = 'eps';
% fig2c  = graph_ci1(time,phi2dot_smp_avg,...
%                    phi2dot_upp,...
%                    phi2dot_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2c  = graph_ci1N(time,MC_phi2dot(samp_pts,:),...
                    phi2dot_upp,...
                    phi2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2cc = graph_type2(time,phi2dot_smp_avg,...
                     time,phi2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2c);
% -----------------------------------------------------------


% plot tower horizontal velocity confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower horizontal velocity';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' velocity (m/s)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_x2_velo'];
gname1 = [num2str(case_name),'__st_x2_velo'];
flag   = 'eps';
% fig2d  = graph_ci1(time,x2dot_smp_avg,...
%                    x2dot_upp,...
%                    x2dot_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2d  = graph_ci1N(time,MC_x2dot(samp_pts,:),...
                    x2dot_upp,...
                    x2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2dd = graph_type2(time,x2dot_smp_avg,...
                     time,x2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2d);
% -----------------------------------------------------------


% plot tower vertical velocity confidence band/statistics
% -----------------------------------------------------------
%gtitle = ' tower vertical velocity';
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
xlab   = ' time (s)';
ylab   = ' velocity (m/s)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__cb_y2_velo'];
gname1 = [num2str(case_name),'__st_y2_velo'];
flag   = 'eps';
% fig2e  = graph_ci1(time,y2dot_smp_avg,...
%                    y2dot_upp,...
%                    y2dot_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2e  = graph_ci1N(time,MC_y2dot(samp_pts,:),...
                    y2dot_upp,...
                    y2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2ee = graph_type2(time,y2dot_smp_avg,...
                     time,y2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2e);
% -----------------------------------------------------------



% plot MC convergence metric
% -----------------------------------------------------------
%gtitle = ' Monte Carlo convergence';
gtitle = ' ';
xlab   = ' number of MC realizations';
ylab   = ' convergence metric';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
%ymin   = 0.22;
%ymax   = 0.25;
gname  = [num2str(case_name),'__MC_conv'];
flag   = 'eps';
fig4   = graph_type1x((1:Ns),MC_conv,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4);
% -----------------------------------------------------------



% plot trailer displacement time-averaged PDF
% -----------------------------------------------------------
%gtitle = ' trailer vertical displacement';
gtitle = ' ';
xlab   = ' time average of (normalized) displacement';
ylab   = ' probability density function';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 2.5;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_y1'];
flag   = 'eps';
fig6a  = graph_bar_curve1(y1_bins,y1_freq,...
                          y1_supp,y1_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6a);
% -----------------------------------------------------------



% plot trailer rotation time-averaged PDF
% -----------------------------------------------------------
%gtitle = ' trailer angular displacement';
gtitle = ' ';
xlab   = ' time average of (normalized) rotation ';
ylab   = ' probability density function';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 4.0;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_phi1'];
flag   = 'eps';
fig6b  = graph_bar_curve1(phi1_bins,phi1_freq,...
                          phi1_supp,phi1_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6b);
% -----------------------------------------------------------


% plot tower rotation time-averaged PDF
% -----------------------------------------------------------
%gtitle = ' tower angular displacement';
gtitle = ' ';
xlab   = ' time average of (normalized) rotation';
ylab   = ' probability density function';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 5.0;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_phi2'];
flag   = 'eps';
fig6c  = graph_bar_curve1(phi2_bins,phi2_freq,...
                          phi2_supp,phi2_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6c);
% -----------------------------------------------------------


% plot tower horizontal displacement time-averaged PDF
% -----------------------------------------------------------
%gtitle = ' tower horizontal displacement';
gtitle = ' ';
xlab   = ' time average of (normalized) displacement';
ylab   = ' probability density function';
xmin   = -4.0;
xmax   =  4.0;
ymin   = 0.0;
ymax   = 1.0;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_x2'];
flag   = 'eps';
fig6d  = graph_bar_curve1(x2_bins,x2_freq,...
                          x2_supp,x2_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6d);
% -----------------------------------------------------------


% plot tower horizontal displacement PDF at time t
% -----------------------------------------------------------
tt1 = 75;
tt2 = 150;
tt3 = 225;
tt4 = 300;

%gtitle = ' tower horizontal displacement';
gtitle = ' ';
xlab   = ' (normalized) displacement';
ylab   = ' probability density function';
xmin   = -4.0;
xmax   =  4.0;
ymin   = 0.0;
ymax   = 1.0;
%xmin   = 'auto';
%xmax   = 'auto';
%ymin   = 'auto';
%ymax   = 'auto';
gname1  = [num2str(case_name),'__pdf_x2_',num2str(tt1)];
gname2  = [num2str(case_name),'__pdf_x2_',num2str(tt2)];
gname3  = [num2str(case_name),'__pdf_x2_',num2str(tt3)];
gname4  = [num2str(case_name),'__pdf_x2_',num2str(tt4)];
flag   = 'eps';
fig6d1 = graph_bar_curve1(x2_bins_t(:,tt1),x2_freq_t(:,tt1),...
                          x2_supp_t(:,tt1), x2_ksd_t(:,tt1),gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig6d2 = graph_bar_curve1(x2_bins_t(:,tt2),x2_freq_t(:,tt2),...
                          x2_supp_t(:,tt2), x2_ksd_t(:,tt2),gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname2,flag);
fig6d3 = graph_bar_curve1(x2_bins_t(:,tt3),x2_freq_t(:,tt3),...
                          x2_supp_t(:,tt3), x2_ksd_t(:,tt3),gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname3,flag);
fig6d4 = graph_bar_curve1(x2_bins_t(:,tt4),x2_freq_t(:,tt4),...
                          x2_supp_t(:,tt4), x2_ksd_t(:,tt4),gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname4,flag);
%close(fig6d);
% -----------------------------------------------------------


% plot trailer vertical displacement entropy
% -----------------------------------------------------------
%gtitle = ' trailer vertical displacement entropy';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' entropy';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__entropy_y1_disp'];
flag   = 'eps';
%fig7a  = graph_type1(time(Njump:end),Sy1,gtitle,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig7a);
% -----------------------------------------------------------


% plot trailer angular displacement entropy
% -----------------------------------------------------------
%gtitle = ' trailer angular displacement entropy';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' entropy';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__entropy_phi1_disp'];
flag   = 'eps';
%fig7b  = graph_type1(time(Njump:end),Sphi1,gtitle,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig7b);
% -----------------------------------------------------------


% plot tower angular displacement entropy
% -----------------------------------------------------------
%gtitle = ' tower angular displacement entropy';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' entropy';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__entropy_phi2_disp'];
flag   = 'eps';
%fig7c  = graph_type1(time(Njump:end),Sphi2,gtitle,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig7c);
% -----------------------------------------------------------



% probability of large lateral vibration
% -----------------------------------------------------------
%gtitle = ' large lateral vibration';
gtitle = ' ';
xlab   = ' time (s)';
ylab   = ' probability';
xmin   = t0;
xmax   = t1;
ymin   = 0.0;
ymax   = 1.0;
gname  = [num2str(case_name),'__prob_x2'];
flag   = 'eps';
fig8a  = graph_type1(time(2:end),1-P_x2(2:end),gtitle,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig8a);
% -----------------------------------------------------------



