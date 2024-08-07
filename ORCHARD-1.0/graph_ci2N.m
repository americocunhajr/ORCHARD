
% -----------------------------------------------------------------
%  graph_ci2N.m
%
%  This functions plots a graph with two curves and an interval
%  of confidence around then.
%
%  input:
%  x1     - x1 data vector (row vector)
%  y1     - y1 data vector (row vector)
%  x2     - x2 data vector (row vector)
%  y2     - y2 data vector (row vector)
%  yupp   - y  upper bound (row vector)
%  ylow   - y  lower bound (row vector)
%  gtitle - graph title
%  leg1   - legend 1
%  leg2   - legend 2
%  leg3   - legend 3
%  xlab   - x axis label
%  ylab   - y axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Mar 23, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_ci2N(x1,y1,x2,y2,yupp,ylow,gtitle,leg1,leg2,leg3,...
                         xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)

                         % check number of arguments
    if nargin < 17
        error('Too few inputs.')
    elseif nargin > 18
        error('Too many inputs.')
    elseif nargin == 17
        flag = 'none';
    end

    % check arguments
    if length(x1) ~= length(y1)
        error('x1 and y1 vectors must be same length')
    end

    if length(x2) ~= length(y2)
        error('x2 and y2 vectors must be same length')
    end
    
    if length(ylow) ~= length(yupp)
        error('ylow and yupp vectors must be same length')
    end
    
    if length(y1) ~= length(ylow)
        error('y1 and ylow vectors must be same length')
    end
    
    % convert to row vectors so fliplr can work
    if find( size(x1) == max(size(x1)) ) < 2
        x1=x1';
    end
    
    if find( size(y1) == max(size(y1)) ) < 2
        y1=y1';
    end
    
    if find( size(x2) == max(size(x2)) ) < 2
        x2=x2';
    end
    
    if find( size(y2) == max(size(y2)) ) < 2
        y2=y2';
    end
    
    if find( size(ylow) == max(size(ylow)) ) < 2
        ylow = ylow';
    end
    
    if find( size(yupp) == max(size(yupp)) ) < 2
        yupp = yupp';
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = plot(x1,y1,'-','Color','b','DisplayName',leg1);
    hold all
    fh2 = plot(x2,y2,'--','Color','r','DisplayName',leg2);
    fh3 = fill([x1 fliplr(x1)],[yupp fliplr(ylow)],[0.9,0.9,0.9],'DisplayName',leg3);
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
    if ( strcmp(xmin,'auto') || strcmp(xmax,'auto') )
        xlim('auto');
    else
        xlim([xmin xmax]);
    end
    
    if ( strcmp(ymin,'auto') || strcmp(ymax,'auto') )
        ylim('auto');
    else
        ylim([ymin ymax]);
    end
    
    set(fh1,'LineWidth',2.5);
    set(fh1,'MarkerSize',2.0);
    set(fh1,'MarkerFaceColor','w');
    set(fh1,'MarkerEdgeColor','k');
    set(fh2,'LineWidth',0.1);
    set(fh2,'MarkerSize',2.0);
    set(fh2,'MarkerFaceColor','w');
    set(fh2,'MarkerEdgeColor','k');
    set(fh3,'LineWidth',1.0);
    set(fh3,'EdgeColor',[0.75,0.75,0.75]);
    set(fh3,'FaceColor',[0.75,0.75,0.75]);
    %leg = legend(leg2,leg1,leg3);
	fh22 = plot(x2,y2(1,:),'--','Color','r','DisplayName',leg1);
    leg = legend([fh1; fh22(1); fh3]);
    set(leg,'FontSize',12);
    %set(leg,'interpreter', 'latex');
    labX = xlabel(xlab,'FontSize',18,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',18,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    uistack(fh3,'top');
    uistack(fh1,'top');
    uistack(fh2,'top');
    
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------
