
% -----------------------------------------------------------------
%  post_orchard_ivp_kl.m
%
%  This script performs the post processing of simulation data
%  that comes from orchard nonlinear dynamics.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 16, 2017
% -----------------------------------------------------------------



% plot trailer vertical displacement PSD
% ...........................................................
gtitle = ' ';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = 2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_y1_disp'];
flag   = 'eps';

xline = [1e0 1e2];
yline = [-20 -60];

 inc = -2;
xinc =  10;
yinc = -33;

%fig3a  = graph_type2_psd(freq,10*log10(psd_Qy1),...
%                         xline,yline,inc,xinc,yinc,...
%                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig3a  = graph_type1(freq,10*log10(psd_Qy1),...
%                     gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3a);
% ...........................................................


% plot trailer angular displacement PSD
% ...........................................................
gtitle = ' ';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = 2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_phi1_disp'];
flag   = 'eps';

xline = [1e0 1e2];
yline = [-20 -60];

 inc = -2;
xinc =  10;
yinc = -33;

%fig3b  = graph_type2_psd(freq,10*log10(psd_Qphi1),...
%                         xline,yline,inc,xinc,yinc,...
%                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig3b  = graph_type1(freq,10*log10(psd_Qphi1),...
%                     gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3b);
% ...........................................................


% plot tower angular displacement PSD
% ...........................................................
gtitle = ' ';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = 2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_phi2_disp'];
flag   = 'eps';

xline = [1e0 1e2];
yline = [-20 -60];

 inc = -2;
xinc =  10;
yinc = -33;


%fig3c  = graph_type2_psd(freq,10*log10(psd_Qphi2),...
%                         xline,yline,inc,xinc,yinc,...
%                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig3c  = graph_type1(freq,10*log10(psd_Qphi2),...
%                     gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3c);
% ...........................................................


% plot tower horizontal displacement PSD
% ...........................................................
gtitle = ' ';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 1.0e-2;
xmax   = 1.0e2;
ymin   = -100;
ymax   = 0;

xline = [1e0 1e2];
yline = [-20 -60];

 inc = -2;
xinc =  10;
yinc = -33;

gname  = [num2str(case_name),'__psd_x2_disp_log'];
flag   = 'eps';
fig3d  = graph_type2_psd(freq,10*log10(psd_Qx2),...
                         xline,yline,inc,xinc,yinc,...
                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
xmin   = 0;
xmax   = 1;
ymin   = -100;
ymax   = 0;
gname  = [num2str(case_name),'__psd_x2_disp'];
fig3dd  = graph_type1(freq,10*log10(psd_Qx2),...
                     gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3d);

line([omega1 omega1],[ymin ymax],'Color','k','LineStyle','--','LineWidth',1.1);
% ...........................................................


% plot tower vertical displacement PSD
% ...........................................................
gtitle = ' ';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_y2_disp'];
flag   = 'eps';

xline = [1e0 1e2];
yline = [-20 -60];

 inc = -2;
xinc =  10;
yinc = -33;

%fig3e  = graph_type2_psd(freq,10*log10(psd_Qy2),...
%                         xline,yline,inc,xinc,yinc,...
%                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig3e  = graph_type1(freq,10*log10(psd_Qy2),...
%                     gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3e);
% ...........................................................

