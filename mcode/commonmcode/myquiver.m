function varargout=myquiver(x,y,vx,vy,varargin)
    sz=10;
    if (isempty(findStringInCell(varargin,'color')))
        quiverColor=[0.5,0.5,0.5];
    else
        quiverColor=varargin{findStringInCell(varargin,'color')+1};
    end
    if (isempty(findStringInCell(varargin,'lw')))        
        lw=1.5;
    else
        lw=varargin{findStringInCell(varargin,'lw')+1};
    end
    
    for n=1:length(x)
        plot([x(n),x(n)+vx(n)],[y(n),y(n)+vy(n)],'Color',quiverColor,'LineWidth',lw);    hold on;
        if (isempty(findStringInCell(varargin,'noRedDots')))
            plot(x(n)+vx(n),y(n)+vy(n),'r.');
        end
        if (abs(vy(n))<1e-8 & abs(vx(n))<1e-8)
            continue;
        end
        if (abs(vy(n))<1e-8)
            if (vx(n)<=0)
                plot([x(n)+vx(n),x(n)+vx(n)+sz],[y(n)+vy(n),y(n)+vy(n)+sz],'Color',quiverColor,'LineWidth',lw);
                plot([x(n)+vx(n),x(n)+vx(n)+sz],[y(n)+vy(n),y(n)+vy(n)-sz],'Color',quiverColor,'LineWidth',lw);
            else
                plot([x(n)+vx(n)-sz,x(n)+vx(n)],[y(n)+vy(n)+sz,y(n)+vy(n)],'Color',quiverColor,'LineWidth',lw);
                plot([x(n)+vx(n)-sz,x(n)+vx(n)],[y(n)+vy(n)-sz,y(n)+vy(n)],'Color',quiverColor,'LineWidth',lw);
            end
        elseif (abs(vx(n))<1e-8)
            if (vy(n)<=0)
                plot([x(n)+vx(n),x(n)+vx(n)+sz*0.5],[y(n)+vy(n),y(n)+vy(n)+sz],'Color',quiverColor,'LineWidth',lw);
                plot([x(n)+vx(n)-sz*0.5,x(n)+vx(n)],[y(n)+vy(n)+sz,y(n)+vy(n)],'Color',quiverColor,'LineWidth',lw);
            else
                plot([x(n)+vx(n),x(n)+vx(n)+sz*0.5],[y(n)+vy(n),y(n)+vy(n)-sz],'Color',quiverColor,'LineWidth',lw);
                plot([x(n)+vx(n)-sz*0.5,x(n)+vx(n)],[y(n)+vy(n)-sz,y(n)+vy(n)],'Color',quiverColor,'LineWidth',lw);
            end            
        end
    end
return