
% -----------------------------------------------------------------
%  plot_orchard_animation.m
%
%  This functions plots an animation for the nonlinear dynamics
%  of an orchard sprayer tower.
%
%  input:
%  time        - (1 x Ndt) time vector (s)
%  y1          - (1 x Ndt) trailer vertical displacement (m)
%  phi1        - (1 x Ndt) trailer anglular displacement (rad)
%  phi2        - (1 x Ndt) tower anglular displacement (rad)
%  B1          - trailer base left arm length (m)
%  B2          - trailer base right arm length (m)
%  L1          - trailer vertical arm length (m)
%  L1          - tower arm length (m)
%  Dtire       - tires diamenter (m)
%  vtrans_km_h - velocity of translation (km/h)
%  vtitle      - video title
%  legend      - legend text
%  xmin        - x axis minimum value
%  xmax        - x axis maximum value
%  ymin        - y axis minimum value
%  ymax        - y axis maximum value
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jun 21, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = post_orchard_animation(time,y1,phi1,phi2,...
                                        B1,B2,L1,L2,Dtire,vtrans_km_h,...
                                        vtitle,vname,legend,xmin,xmax,ymin,ymax)
    
    % check number of arguments
    if nargin < 17
        error('Too few inputs.')
    elseif nargin > 18
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(y1)
        error('vectors time and y1 must have same length')
    end
    
    if length(time) ~= length(phi1)
        error('vectors time and phi1 must have same length')
    end
    
    if length(time) ~= length(phi2)
        error('vectors time and phi2 must have same length')
    end
    
    if B1 <= 0.0
        error('B1 must be positive')
    end
    
    if B2 <= 0.0
        error('B2 must be positive')
    end
    
    if L1 <= 0.0
        error('L1 must be positive')
    end
    
    if L2 <= 0.0
        error('L2 must be positive')
    end
    
    if Dtire <= 0.0
        error('Dtire must be positive')
    end
    
    % convert to row vector (if necessary)
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(y1) == max(size(y1)) ) < 2
        y1=y1';
    end
    
    if find( size(phi1) == max(size(phi1)) ) < 2
        phi1=phi1';
    end
    
    if find( size(phi2) == max(size(phi2)) ) < 2
        phi2=phi2';
    end
    
    
    % number of time steps
    Ndt = length(time);
    
    % coordinates of the trailer center of mass
    xcm = zeros(1,Ndt);
    ycm = y1;

    % coordinates of the trailer border
    xtrailer_l = xcm - B1*cos(phi1);
    ytrailer_l = ycm - B1*sin(phi1);
    xtrailer_r = xcm + B2*cos(phi1);
    ytrailer_r = ycm + B2*sin(phi1);

    % coordinates of the pivot
    xpivot = xcm - L1*sin(phi1);
    ypivot = ycm + L1*cos(phi1);

    % tower center of mass coordinates
    xtower = xpivot - L2*sin(phi2);
    ytower = ypivot + L2*cos(phi2);

    % coordinates of the time counter
    xtime = 0.4*xmax;
    ytime = 0.9*ymax;
    
    % coordinates of velocity indicator
    xvelo = 0.8*xmin;
    yvelo = 0.9*ymax;
    
    % coordinates of legend
    xleg = 0.8*xmin;
    yleg = 0.8*ymax;
    
    % soil heigth (m)
    ysoil = 0.25*(ymax-ymin);

    
    % loop to construct the video
    for n=1:Ndt
        
        % initialize video frame
        fig = figure(1000);
        
        % define frame properties
        set(gcf,'color','white');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'zColor',[.3 .3 .3]);
        set(gca,'FontName','Helvetica');
        set(gca,'FontSize',16);
        set(gca,'Box','on');
        
        
        % draw soil border
        fh01 = plot([xmin xmax],[(ymin+ysoil) (ymin+ysoil)],...
                                               '-k','LineWidth',2);
        
        hold on
        
        % draw soil
        fh02  = rectangle('Position',[xmin,ymin,(xmax-xmin),ysoil],...
                         'FaceColor',[0.95,0.95,0.95]);

        % draw trailer base right arm
        fh03 = plot([xcm(n) xtrailer_r(n)],[ycm(n) ytrailer_r(n)],...
                                                   '-r','LineWidth',2);
        
        % draw trailer base left arm
        fh04 = plot([xcm(n) xtrailer_l(n)],[ycm(n) ytrailer_l(n)],...
                                                   '-r','LineWidth',2);

        % draw trailer vertical arm
        fh05 = plot([xcm(n) xpivot(n)],[ycm(n) ypivot(n)],...
                                                   '-r','LineWidth',2);

        % draw tower
        fh06 = plot([xpivot(n) xtower(n)],[ypivot(n) ytower(n)],...
                                                   '-b','LineWidth',6);

        % draw left wheel
        fh07 = rectangle('Position',[xtrailer_l(n)-0.1*Dtire,...
                                     ytrailer_l(n)-0.5*Dtire,...
                                     0.1*Dtire,Dtire],...
                         'FaceColor','k');

        % draw right wheel
        fh08 = rectangle('Position',[xtrailer_r(n),...
                                     ytrailer_r(n)-0.5*Dtire,...
                                     0.1*Dtire,Dtire],...
                         'FaceColor','k');
        
        % define video frame limits
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        
        % display time counter
        text(xtime,ytime,['time = ',num2str(time(n),'%.1f'),' sec'],...
                          'Color','k',...
                          'FontName','Helvetica',...
                          'FontSize',14);
                      
        % display velocity indicator
        text(xvelo,yvelo,['velocity = ',num2str(vtrans_km_h,'%.1f'),' km/h'],...
                          'Color','k',...
                          'FontName','Helvetica',...
                          'FontSize',14);
        
        % display legend text
        text(xleg,yleg,legend,...
                       'Color','k',...
                       'FontName','Helvetica',...
                       'FontSize',14);

                      
        % define video title
        title(vtitle,'FontSize',16,'FontName','Helvetica');
        
        hold off
        
        % pause exibition
        if time(n) == 0.0
            pause
        end
        
        if abs(time(n) - 3.0) < 0.01
            pause
        end
        
        if abs(time(n) - 7.0) < 0.01
            pause
        end
        
        if abs(time(n) - 11.0) < 0.01
            pause
        end
        
        if abs(time(n) - 15.0) < 0.01
            pause
        end
        
        if abs(time(n) - 18.0) < 0.01
            pause
        end
        
        if abs(time(n) - 22.0) < 0.01
            pause
        end
        
        if abs(time(n) - 26.0) < 0.01
            pause
        end
        
        if abs(time(n) - 30.0) < 0.01
            pause
        end

    end
    
end
% -----------------------------------------------------------------
