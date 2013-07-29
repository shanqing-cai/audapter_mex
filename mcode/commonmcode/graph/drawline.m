function drawline(x,varargin)
% Drawline(x,varargin)
% VARGIN: 
%   'dir','vert'        ('horz' gives horizontal lines) 
%   'color',[0 0 0]     (default black)
%   'lim',[min max]     (default: existing axis limits) 
color=[0 0 0];
lim=[]; 
dir='vert'; 

[r,c]=size(x);
if(r>c)
   x=x';
end;

vararginoptions(varargin,{'dir','color','lim'});

if (isempty(lim))
    if (strcmp(dir,'vert'))
        lim=get(gca,'YLim');
    else 
        lim=get(gca,'XLim');
    end;
end;

x1=[x;x];
y=[ones(1,length(x));ones(1,length(x))];
y(1,:)=y(1,:).*lim(1);
y(2,:)=y(2,:).*lim(2);

if (strcmp(dir,'vert')) 
    line(x1,y,'Color',color);
else 
    line(y,x1,'Color',color);
end;
    
