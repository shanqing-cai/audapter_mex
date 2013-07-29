function varargout=scatterplot(x,y,varargin)
%  function scatterplot(x,y,varargin)
% Description
%  varargin:
%   Format options (for all symbols) 
%       'markertype',{o s v ^...}   
%       'markercolor',[r g b] 
%       'markerfill',[r g b]
%       'markersize',size
%   Other commands:
%       'CAT',cat        : Stucture with fields with one entry per category (split by) (all fields are 
%                 cell arrays with length of the number of categories 
%       'subset',indx    : plots only a subset of the data
%       'split',splitby  : splits data by vaiable splitby
%       'color',x      : Uses current color map to shade the circle in
%                           color corresponding to the range.
%       'colormap', CM: sets colormap
%       'leg',{}/'auto'  : specified legend or auto
%       'regression','linear':       Puts a linear regression line per category in 
%       'corr':         Prints correlation onto the graph 
%       'draworig':       Draws origin into the graph
% Joern Diedrichsen (jdiedric@bme.jhu.edu)
% v.1.1 9/18/05
% v.1.2 12/14/05: Support for data sets without variation is added 

[Nx n] = size(y);
if (n>1)  
    error('data has to be a row-vector');
end

% Set defaults for all plots 
markertype = 'o';
markercolor=[0 0 0];
markerfill=[1 1 1];
colormap=hot;
markersize=4;
CAT.markertype= {'o','o','o','o','o','o','s','s','s','s','s','s','v','v','v','v','v','v'};
CAT.markercolor={[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1]};
CAT.markerfill={[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 0 1],[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 0 1],[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 0 1]};
draworig=0;
corr=0;
regression=[];
leg=[];
leglocation='SouthEast';
% Deal with the varargin's 
split=[];subset=ones(size(y,1),1);color=[];
variables={'markertype','markercolor','markerfill','markersize','CAT','subset','split','leg','leglocation','color','regression','colormap'};
flags={'draworig','corr'};

vararginoptions(varargin,variables,flags);
F.markertype = markertype;
F.markercolor=markercolor;
F.markerfill=markerfill;
F.markersize=markersize;
F.colormap=colormap;

if (ischar(CAT))
    switch (CAT)
        case '4x2x2'
            CAT.markertype= {'o','o','o','o','o','o','o','o','s','s','s','s','s','s','s','s'};
            CAT.markercolor={[0 0 0],[1 0 0],[0 1 0],[0 0 1],[0 0 0],[1 0 0],[0 1 0],[0 0 1],[0 0 0],[1 0 0],[0 1 0],[0 0 1],[0 0 0],[1 0 0],[0 1 0],[0 0 1]};
            CAT.markerfill={[1 1 1],[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 1 0],[0 0 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 1 0],[0 0 1]};
        case '3x2x2'
            CAT.markertype= {'o','o','o','o','o','o','s','s','s','s','s','s'};
            CAT.markercolor={[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1],[0 0 0],[1 0 0],[0 0 1]};
            CAT.markerfill={[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 0 1],[1 1 1],[1 1 1],[1 1 1],[0 0 0],[1 0 0],[0 0 1]};
        case '2x8'
            CAT.markertype= {'o','o','+','+','*','*','x','x','s','s','<','<','>','>','^','^'};
            CAT.markercolor={[1 0 0],[0 1 0],[1 0 0],[0 1 0],[1 0 0],[0 1 0],[1 0 0],[0 1 0],...
                    [1 0 0],[0 1 0],[1 0 0],[0 1 0],[1 0 0],[0 1 0],[1 0 0],[0 1 0]};
        otherwise
            error(['CAT: Undefined Map' CAT]);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-code cell array  
if ~isempty(split)
    [split,split_conv]=fac2int(split);
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply goodindx-selection for subset
goodindx=find(subset);
y=y(goodindx,:);
x=x(goodindx,:);
if (~isempty(split))
    split=split(goodindx,:);
end;
if (~isempty(color))
    color=color(goodindx,:);
    colorz=round((size(colormap,1)-1)*(color-min(color))./(max(color)-min(color))+1);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate scale 
xmin = min(x);
xmax = max(x);
dx = (xmax-xmin)/20;
if (dx==0) dx=xmax/10;end;
xlims = [(xmin-dx) (xmax+dx)];

ymin = min(y);
ymax = max(y);
dy = (ymax-ymin)/20;
if (dy==0) dy=ymax/10;end;
ylims = [(ymin-dy) (ymax+dy)];


repl_state=get(gca,'NextPlot');
if (strcmp(repl_state,'replace'))
    cla;
end;
set(gca,'NextPlot','add');
set(gca,'Box','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with split 
if (~isempty(split))
    D=pidata(split,[x y]);
    [numsplitcat,dummy]=size(D);
    for c=1:numsplitcat
        if (~isempty(CAT))
            fields=fieldnames(CAT);
            for f=1:length(fields)
                fiel=getfield(CAT,fields{f});
                %F.(fields{f})=fiel{spcat};
                F=setfield(F,fields{f},fiel{c});
            end;
        end;
        h(c)=plot(D{c,2}(:,1),D{c,2}(:,2),'ko');
        set(h(c),'Marker',F.markertype,'MarkerSize',F.markersize,...
                'MarkerEdgeColor',F.markercolor,'MarkerFaceColor',F.markerfill);
        if (~isempty(regression))
            beta(:,c)=doregress(D{c,2}(:,1),D{c,2}(:,2),F.markercolor);
            varargout={beta};
        end;
    end; 
else
    if (isempty(color))
        set(plot(x,y,'k.'),'Marker',F.markertype,'MarkerSize',F.markersize,...
            'MarkerEdgeColor',F.markercolor,'MarkerFaceColor',F.markerfill);
    else
        for i=1:length(x)
            set(plot(x(i),y(i),'k.'),'Marker',F.markertype,'MarkerSize',F.markersize,...
                'MarkerEdgeColor',colormap(colorz(i),:),'MarkerFaceColor',colormap(colorz(i),:));hold on;
            end;
        hold off;
    end;
    if (~isempty(regression))
        beta=doregress(x,y,[0 0 0]);
        varargout={beta};
    end;
end;
axis([xlims ylims]);
if (draworig)
    draworigin;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do legend 
% Add legend if necessary 
if (~isempty(split))
    Split_groups=vertcat(D{:,1});
    plotlegend(h,leg,Split_groups,split_conv,leglocation);
end;
set(gca,'NextPlot',repl_state);

function b=doregress(x,y,color);
X=[ones(length(x),1) x];
b=inv(X'*X)*X'*y;
Xa=[min(x) max(x)];
Ya=b(2)*Xa+b(1);
h=line(Xa,Ya);
set(h,'Color',color);