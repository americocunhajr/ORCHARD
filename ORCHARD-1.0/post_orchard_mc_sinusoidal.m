
% -----------------------------------------------------------------
%  post_orchard_mc_sinusoidal.m
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
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig1a  = graph_ci1N(time,y1(samp_pts,:),...
                    y1_upp,...
                    y1_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
                
fig1aa = graph_type2(time,y1_smp_avg,...
                     time,y1_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig1aaa = graph_type4(time,y1_smp_avg,...
                      time,y1_std,...
                      time,y1_skew,...
                      time,y1_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1a);
% -----------------------------------------------------------


% plot trailer angular displacement confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig1b  = graph_ci1N(time,phi1(samp_pts,:),...
                    phi1_upp,...
                    phi1_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1bb = graph_type2(time,phi1_smp_avg,...
                     time,phi1_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig1bbb = graph_type4(time,phi1_smp_avg,...
                      time,phi1_std,...
                      time,phi1_skew,...
                      time,phi1_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1b);
% -----------------------------------------------------------



% plot tower angular displacement confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig1c  = graph_ci1N(time,phi2(samp_pts,:),...
                     phi2_upp,...
                     phi2_low,...
                     gtitle,leg1,leg2,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1cc = graph_type2(time,phi2_smp_avg,...
                     time,phi2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig1ccc = graph_type4(time,phi2_smp_avg,...
                      time,phi2_std,...
                      time,phi2_skew,...
                      time,phi2_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1c);
% -----------------------------------------------------------



% plot tower horizontal displacement confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
xlab   = ' time (s)';
ylab   = ' displacement (m)';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
%ymin   = -0.2;
%ymax   = 0.2;
ymin1   = 'auto';
ymax1   = 'auto';
gname  = [num2str(case_name),'__cb_x2_disp'];
gname1 = [num2str(case_name),'__st_x2_disp'];
flag   = 'eps';
% fig1d  = graph_ci1(time,x2_smp_avg,...
%                    x2_upp,...
%                    x2_low,...
%                    gtitle,leg2,leg3,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1d  = graph_ci1N(time,x2(samp_pts,:),...
                    x2_upp,...
                    x2_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1dd = graph_type2(time,x2_smp_avg,...
                     time,x2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin1,ymax1,gname1,flag);
%fig1ddd = graph_type4(time,x2_smp_avg,...
%                      time,x2_std,...
%                      time,x2_skew,...
%                      time,x2_kurt,...
%                      gtitle,leg3,leg4,leg5,leg6,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1d);
% -----------------------------------------------------------



% plot tower vertical displacement confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig1e  = graph_ci1N(time,y2(samp_pts,:),...
                    y2_upp,...
                    y2_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1ee = graph_type2(time,y2_smp_avg,...
                     time,y2_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig1eee = graph_type4(time,y2_smp_avg,...
                      time,y2_std,...
                      time,y2_skew,...
                      time,y2_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig1e);
% -----------------------------------------------------------


% plot trailer vertical velocity confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig2a  = graph_ci1N(time,y1dot(samp_pts,:),...
                    y1dot_upp,...
                    y1dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2aa = graph_type2(time,y1dot_smp_avg,...
                     time,y1dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig2aaa = graph_type4(time,y1dot_smp_avg,...
                      time,y1dot_std,...
                      time,y1dot_skew,...
                      time,y1dot_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2a);
% -----------------------------------------------------------



% plot trailer angular velocity confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig2b  = graph_ci1N(time,Qphi1dot(samp_pts,:),...
                    phi1dot_upp,...
                    phi1dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2bb = graph_type2(time,phi1dot_smp_avg,...
                     time,phi1dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig2bbb = graph_type4(time,phi1dot_smp_avg,...
                      time,phi1dot_std,...
                      time,phi1dot_skew,...
                      time,phi1dot_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2b);
% -----------------------------------------------------------



% plot tower angular velocity confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig2c  = graph_ci1N(time,Qphi2dot(samp_pts,:),...
                    phi2dot_upp,...
                    phi2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2cc = graph_type2(time,phi2dot_smp_avg,...
                     time,phi2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig2ccc = graph_type4(time,phi2dot_smp_avg,...
                      time,phi2dot_std,...
                      time,phi2dot_skew,...
                      time,phi2dot_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2c);
% -----------------------------------------------------------


% plot tower horizontal velocity confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig2d  = graph_ci1N(time,Qx2dot(samp_pts,:),...
                    x2dot_upp,...
                    x2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2dd = graph_type2(time,x2dot_smp_avg,...
                     time,x2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig2ddd = graph_type4(time,x2dot_smp_avg,...
                      time,x2dot_std,...
                      time,x2dot_skew,...
                      time,x2dot_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2d);
% -----------------------------------------------------------


% plot tower vertical velocity confidence band/statistics
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' realizations';
leg2   = [' ',num2str(Pc),'% prob.'];
leg3   = ' mean';
leg4   = ' std. dev.';
leg5   = ' skewness';
leg6   = ' kurtosis';
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
fig2e  = graph_ci1N(time,Qy2dot(samp_pts,:),...
                    y2dot_upp,...
                    y2dot_low,...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2ee = graph_type2(time,y2dot_smp_avg,...
                     time,y2dot_std,...
                     gtitle,leg3,leg4,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
fig2eee = graph_type4(time,y2dot_smp_avg,...
                      time,y2dot_std,...
                      time,y2dot_skew,...
                      time,y2dot_kurt,...
                      gtitle,leg3,leg4,leg5,leg6,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname1,flag);
%close(fig2e);
% -----------------------------------------------------------



% plot MC convergence metric
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' number of MC realizations';
ylab   = ' convergence metric';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__MC_conv'];
flag   = 'eps';
fig4   = graph_type1x((1:Ns),MC_conv,gtitle,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4);
% -----------------------------------------------------------



% plot trailer phase portrait
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' nominal';
leg2   = ' mean value';
xlab   = ' displacement (m)';
ylab   = ' velocity (m/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_y1'];
flag   = 'eps';
%fig5a  = graph_type2(Qy1_nominal(Nss:end),Qy1dot_nominal(Nss:end),...
%                     y1_smp_avg(Nss:end),y1dot_smp_avg(Nss:end),...
%                     gtitle,leg1,leg2,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5a);
% -----------------------------------------------------------


% plot trailer angular phase portrait
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' nominal';
leg2   = ' mean value';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi1'];
flag   = 'eps';
%fig5b  = graph_type2(Qphi1_nominal(Nss:end),Qphi1_dot_nominal(Nss:end),...
%                     phi1_smp_avg(Nss:end),phi1_dot_smp_avg(Nss:end),...
%                     gtitle,leg1,leg2,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5b);
% -----------------------------------------------------------


% plot tower angular phase portrait
% -----------------------------------------------------------
gtitle = ' ';
leg1   = ' nominal';
leg2   = ' mean value';
xlab   = ' rotation (rad)';
ylab   = ' angular velocity (rad/s)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__phase_port_phi2'];
flag   = 'eps';
%fig5c  = graph_type2(Qphi2_nominal(Nss:end),Qphi2_dot_nominal(Nss:end),...
%                     phi2_smp_avg(Nss:end),phi2_dot_smp_avg(Nss:end),...
%                     gtitle,leg1,leg2,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig5c);
% -----------------------------------------------------------




% plot trailer displacement (steady state) time average PDF
% -----------------------------------------------------------
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



% plot trailer rotation (steady state) time average PDF
% -----------------------------------------------------------
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


% plot tower rotation (steady state) time average PDF
% -----------------------------------------------------------
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


% plot tower horizontal displacement (steady state) time average PDF
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' time average of (normalized) displacement';
ylab   = ' probability density function';
xmin   = -5.0;
xmax   =  5.0;
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
idx = Ndt;

gtitle = ' ';
xlab   = ' displacement (normalized)';
ylab   = ' probability density function';
xmin   = -4;
xmax   =  4;
ymin   = 0;
ymax   = 1;
gname  = [num2str(case_name),'__pdf_x2_t',num2str(idx)];
flag   = 'eps';
fig66  = graph_bar_curve1leg(x2_bins_t(:,idx),x2_freq_t(:,idx),...
                             x2_supp_t(:,idx),x2_ksd_t(:,idx),...
                             time(idx),gtitle,...
                             xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig66);
% -----------------------------------------------------------


% 
% % plot trailer vertical displacement entropy
% % -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' time (s)';
% ylab   = ' entropy';
% xmin   = t0;
% xmax   = t1;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__entropy_y1_disp'];
% flag   = 'eps';
% fig7a  = graph_type1(time(Njump:end),Sy1,gtitle,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
% %close(fig7a);
% % -----------------------------------------------------------
% 
% 
% % plot trailer angular displacement entropy
% % -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' time (s)';
% ylab   = ' entropy';
% xmin   = t0;
% xmax   = t1;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__entropy_phi1_disp'];
% flag   = 'eps';
% fig7b  = graph_type1(time(Njump:end),Sphi1,gtitle,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
% %close(fig7b);
% % -----------------------------------------------------------
% 
% 
% % plot tower angular displacement entropy
% % -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' time (s)';
% ylab   = ' entropy';
% xmin   = t0;
% xmax   = t1;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__entropy_phi2_disp'];
% flag   = 'eps';
% fig7c  = graph_type1(time(Njump:end),Sphi2,gtitle,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
% %close(fig7c);
% % -----------------------------------------------------------



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
fig8a  = graph_type1(time(2:end),1-prob_x2(2:end),gtitle,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig8a);
% -----------------------------------------------------------


