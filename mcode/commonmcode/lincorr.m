function [k,r2,p]=lincorr(x,y,varargin)
    if (size(x,2)>size(x,1))
        x=x';
    end
    if (size(y,2)>size(y,1));
        y=y';
    end
    
    [b,bint,r,rint,stats]=regress(y,[ones(size(x)),x]);
    k=b;
    r2=stats(1);
    p=stats(3);
    
    if (~isempty(findStringInCell(varargin,'plot')))
        xs=varargin{findStringInCell(varargin,'plot')+1};
        ys=varargin{findStringInCell(varargin,'plot')+2};
        xoffset=varargin{findStringInCell(varargin,'plot')+3};
        yoffset=varargin{findStringInCell(varargin,'plot')+4};
        color=varargin{findStringInCell(varargin,'plot')+5};
    
        plot(xs,xs*k(2)+k(1),'Color',color);
        text(xs(1)+xoffset*range(xs),ys(1)+yoffset*range(ys),['r^2=',num2str(r2),'; p=',num2str(p)],'Color',color);
    end
    
return