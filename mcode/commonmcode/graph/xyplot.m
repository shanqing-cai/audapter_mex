function [X,Y,XERR,YERR]=xyplot(xvar,yvar,cat,varargin)
% Synopsis
%  [X,Y,Xerr,Yerr]=lineplot(xvar,yvar,cat,varargin) 
% Description
% plots a variable against another using x and y error bars 
% data series is categoriezed (same symbols - line) or 
% split (different symbols)
%  X: variable 1 [N*1]
%  Y: varaible 2 [N*1]
%
%  varargin:
%   Format options (for all lines) 
%       'markersize',size    : Size of the marker symbols 
%       'markertype',type    : Type ('o','s','v','^') of the marker 
%       'markercolor',[r g b]: Color of the marker edge 
%       'markerfill',[r g b] : Color of the marker fill
%       'linecolor',[r g b]  : Color of the Line  
%       'linewidth',width    : Width of the Line 
%       'CAT', CAT           : Structure with fields with one entry per category (split by) 
%                               Field can be all of the above (see online
%                               doc)
%    Predetermined styles 
%       'style_symbols4*2': Square, circle,Triangle (up/down), in white and
%                     black
%       'style_thickline':  Thick lines in different colors, square symbols, and error
%                     bars,different colors 
%   Data processing options 
%       'plotfcn_x'   : function over data of what should be plotted: default 'mean'
%       'plotfcn_y'   : function over data of what should be plotted: default 'mean'
%       'errorfcn_x'  : function over data to determine size of error bars: default 'stderr'
%       'errorfcn_y'  : function over data to determine size of error bars: default 'stderr'
%       'split',var   : Variable to split the data by. Seperate lines are
%                        drawn per value of the split-var 
%       'subset'      : Plots only a subset of the data
%       'leg'         : Legend, either cell with names or 'auto'
%       'leglocation','north' : Legend location          
% v1.0: 8/7/07: Joern diedrichsen 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set defaults for all plots 
F.linecolor=[0.8 0.8 0.8];
F.linewidth=1;
F.linestyle='-';
F.fillcolor=[0.8 0.8 0.8];
F.markertype = 'o';
F.markercolor=[0 0 0];
F.markerfill=[0 0 0];
F.markersize=4;
F.flip=0;
F.errorbars='plusminus';
F.errorwidth=1;
F.errorcolor=[0.8 0.8 0.8];
gap=[1 0.7 0.5 0.5];
leg=[];
leglocation='SouthEast';
CAT.markertype={'o','o','^','^','s','s','v','v'};
CAT.markerfill={[0 0 0],[1 1 1],[0 0 0],[1 1 1],[0 0 0],[1 1 1],[0 0 0],[1 1 1]};

plotfcn_x='mean';
plotfcn_y='mean';
errorfcn_x='stderr';
errorfcn_y='stderr';

numxvars=size(xvar,2);
XTickLabel=0;
XCoord='auto';
split=[];
numsplitvars=0;
subset=ones(size(xvar,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with the varargin's 
c=1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'XTickLabel','XCoord','plotfcn','errorfcn','errorfcn_up','errorfcn_down','CAT','leg','leglocation','split','subset','transformfcn','plotfcn_x','plotfcn_y'}
            eval([varargin{c} '=varargin{c+1};']);
            c=c+2;
        case {'markertype','markerfill','linecolor','linewidth','markercolor','markersize','errorcolor','errorwidth'}
            eval(['F.' varargin{c} '=varargin{c+1};']);
            c=c+2;
        case 'avrgcorr'
            plotfcn='mean(fisherz(x))';
            errorfcn='stderr(fisherz(x))';
            transformfcn='fisherinv(x)';
            type='avrgcorr';c=c+1;
        case 'style_symbols4*2'
            CAT.markertype={'o','o','^','^','s','s','v','v'};
            CAT.markerfill={[0 0 0],[1 1 1],[0 0 0],[1 1 1],[0 0 0],[1 1 1],[0 0 0],[1 1 1]};c=c+1;
        case 'style_thickline'
            CAT.linecolor={[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0]};
            CAT.markercolor={[0 0 1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0]};
            CAT.errorcolor={[0 0  1],[1 0 0],[0 1 0],[1 0 1],[0 1 1],[0.7 0.7 0.7],[1 1 0]};
            F.markertype='s';
            F.markersize=2.5;
            F.linewidth=2.5;
            F.errorwidth=1;
            c=c+1;
        case 'style_thickline2x3'
            CAT.linecolor={[0 0 1],[0 0 1],[1 0 0],[1 0 0],[0 1 0],[0 1 0]};
            CAT.markercolor={[0 0 1],[0 0 1],[1 0 0],[1 0 0],[0 1 0],[0 1 0]};
            CAT.errorcolor={[0 0 1],[0 0 1],[1 0 0],[1 0 0],[0 1 0],[0 1 0]};
            CAT.linestyle={'-','--','-','--','-','--'};
            F.markertype='s';
            F.markersize=2.5;
            F.linewidth=2.5;
            F.errorwidth=1;
            c=c+1;
        case 'style_thickline3x3'
            CAT.linecolor={[0 0 1],[0 0 1],[0 0 1],[1 0 0],[1 0 0],[1 0 0],[0 1 0],[0 1 0],[0 1 0]};
            CAT.markercolor={[0 0 1],[0 0 1],[0 0 1],[1 0 0],[1 0 0],[1 0 0],[0 1 0],[0 1 0],[0 1 0]};
            CAT.errorcolor={[0 0 1],[0 0 1],[0 0 1],[1 0 0],[1 0 0],[1 0 0],[0 1 0],[0 1 0],[0 1 0]};
            CAT.linestyle={'-','--',':','-','--',':','-','--',':'};
            CAT.markertype={'s','v','o','s','v','o','s','v','o'};
            F.markersize=2.5;
            F.linewidth=2.5;
            F.errorwidth=1;
            c=c+1;
        case 'style_shadeline'
            CAT.linecolor={[0.2 0.2 0.2],[0.3 0.3 0.3],[0.4 0.4 0.4],[0.5 0.5 0.5],[0.6 0.6 0.6],[0.7 0.7 0.7],[0.8 0.8 0.8]};
            F.markertype='none';
            F.linewidth=2;
            F.errorwidth=3;
            c=c+1;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recode cat and split variable if they are strings 
if ~isempty(split)
    [split,split_conv]=fac2int(split);
end;    
if ~isempty(cat)
    [cat,cat_conv]=fac2int(cat);
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with selection (subset) variable 
goodindx=find(subset);
if (isempty(goodindx))
    return;
end;
yvar=yvar(goodindx,:);
xvar=xvar(goodindx,:);
if (~isempty(cat))
    cat=cat(goodindx,:);
end;
if (~isempty(split))
    split=split(goodindx,:);
end;

numsplitvars=size(split,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate the mean plot and errorbar size for each category 
[X,R,C]=pivottable(split,cat,xvar,plotfcn_x);
[Y,R,C]=pivottable(split,cat,yvar,plotfcn_y);
if (~isempty(errorfcn_x))
    [Xerr,R,C]=pivottable(split,cat,xvar,errorfcn_x);
end;
if (~isempty(errorfcn_y))
    [Yerr,R,C]=pivottable(split,cat,yvar,errorfcn_y);
end;
numsplitcat=size(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate scale 
xmin = min(X(:)-Xerr(:));
xmax = max(X(:)+Xerr(:));
ymin = min(Y(:)-Yerr(:));
ymax = max(Y(:)+Yerr(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale axis for vertical or horizontal boxes.
% deal with hold on option.
holding=get(gca,'NextPlot');
if (strcmp(holding,'add'))
    xlim_old=get(gca,'XLim');
    ylim_old=get(gca,'YLim');
    xlims=[min(xmin,xlim_old(1)) max(xmax,xlim_old(2))];
    ylims=[min(ymin,ylim_old(1)) max(ymax,ylim_old(2))];
else
    cla
    xlims=[xmin xmax];
    ylims=[ymin ymax];
    set(gca,'NextPlot','add');
end;
set(gca,'Box','off');
axis([xlims ylims]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform plotting 
for c=1:numsplitcat
    if (~isempty(CAT))
        fields=fieldnames(CAT);
        for f=1:length(fields)
            fiel=getfield(CAT,fields{f});
            size_fiel=length(fiel);
            if c>size_fiel
                warning('More categories in split variable than in format structure. Recycling formats!');
            end;
            F=setfield(F,fields{f},fiel{mod(c-1,size_fiel)+1});
        end;
    end;
    h(c)=plot(X(c,:)',Y(c,:)');
    set(h(c),'Color',F.linecolor,'LineWidth',F.linewidth,'LineStyle',F.linestyle,'Marker',F.markertype,'MarkerSize',F.markersize,'MarkerEdgeColor',F.markercolor,'MarkerFaceColor',F.markerfill); 
    if (strcmp(F.errorbars,'plusminus'))
        errorbars(X(c,:)',Y(c,:)',Yerr(c,:)','linecolor',F.errorcolor,'linewidth',F.errorwidth,'orientation','vert','error_dir','both');
        errorbars(X(c,:)',Y(c,:)',Xerr(c,:)','linecolor',F.errorcolor,'linewidth',F.errorwidth,'orientation','horz','error_dir','both');
    % elseif (strcmp(F.errorbars,'ellipse'))
        % To be implemented 
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do legend 
if (~isempty(split))
    plotlegend(h,leg,R,split_conv,leglocation);
end;


