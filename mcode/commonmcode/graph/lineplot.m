function [x_coord,PLOT,ERROR]=lineplot(xvar,y,varargin)
% Synopsis
%  [PLOT,ERROR]=lineplot(xvar,y,varargin)
% Description
%  xvar: independent variables [N*c], with c>1 a hierarchical grouping is used
%  Y: dependent variable [N*1]
%     if Y is a N*p varaible, then different lines are plotted for
%     different variables (like split). In this case, split can't be used anymore
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
%       'plotfcn'   : function over data of what should be plotted: default 'mean'
%       'errorfcn'  : function over data to determine size of error bars: default 'stderr'
%                     if just one error function is given, we assume
%                     symmetric error bars
%       'errorfcn_up': If given, error bars can be asymmetric
%       'errorfcn_down':
%       'transformfcn': All plots-elements (including error bars are
%                       transformed before plotting
%       'split',var   : Variable to split the data by. Seperate lines are
%                        drawn per value of the split-var
%       'subset'      : Plots only a subset of the data
%       'leg'         : Legend, either cell with names or 'auto'
%       'leglocation','north' : Legend location
%  EXAMPLE: for averaging of correltation after fisher-z transformation and
%       concurrently inverse fisher-z trasnform + correct standard errors:
%       lineplot(CORR,X,'plotfcn','mean(fisherz(x))',
%                       'errorfcn','stderr(fisherz(x))'
%                       'transformfcn','fisherinv';
%                         .....
%                       SHORTCUT: 'avrgcorr'
% v1.1: 12/14/05: Support for data with variation (constant data) is added
%       12/22/05: Warning instead of error when number of categories
%       exceeds the number of formats, recycling formats
if (nargin==1 & length(y(:))==1 & ishandle(y)), resizefcn(y); return; end

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
plotfcn='mean';
errorfcn='stderr';
numxvars=size(xvar,2);
XTickLabel=0;
XCoord='auto';
split=[];numsplitvars=0;
goodindx=[1:size(y,1)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with the varargin's
c=1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'gap','XTickLabel','XCoord','plotfcn','errorfcn','errorfcn_up','errorfcn_down','CAT','leg','leglocation'}
            eval([varargin{c} '=varargin{c+1};']);
            c=c+2;
        case {'markertype','markerfill','linecolor','linewidth','markercolor','markersize','errorcolor','errorwidth'}
            eval(['F.' varargin{c} '=varargin{c+1};']);
            c=c+2;
        case 'subset'
            goodindx=find(varargin{c+1});
            c=c+2;
        case 'split'
            split=varargin{c+1};c=c+2;
            [dummy,numsplitvars]=size(split);
        case 'avrgcorr'
            plotfcn='mean(fisherz(x))';
            errorfcn='stderr(fisherz(x))';
            transformfcn='fisherinv(x)';
            type='avrgcorr';c=c+1;
        case 'transformfcn'
            transformfcn=varargin{c+1};c=c+2;
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
% re-code cell array
[xvar,xvar_conv]=fac2int(xvar);
if ~isempty(split)
    [split,split_conv]=fac2int(split);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with selection (subset) variable
y=y(goodindx,:);
xvar=xvar(goodindx,:);
if (~isempty(split))
    split=split(goodindx,:);
end;
if (isempty(y))
    return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with y-variables which have more than one column
[Nx p] = size(y);
if (p>1)
    y=reshape(y,prod(size(y)),1);
    xvar=repmat(xvar,p,1);
    for i=1:p
        split=[split;ones(Nx,1)*i];
    end;
    [split,split_conv]=fac2int(split);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate the mean plot and errorbar size for each category
[PLOT,R,C]=pivottable(split,xvar,y,plotfcn);
if (all(isnan(PLOT(:))))
    error('no data to plot - all is nan');
end;
if (exist('errorfcn_up'))
    [ERROR_UP,R,C]=pivottable(split,xvar,y,errorfcn_up);
end;
if (exist('errorfcn_down'))
    [ERROR_DOWN,R,C]=pivottable(split,xvar,y,errorfcn_down);
end;
if (~isempty(errorfcn) & ~exist('ERROR_UP') & ~ exist('ERROR_DOWN'))
    [ERROR,R,C]=pivottable(split,xvar,y,errorfcn);
    ERROR_UP=ERROR;
    ERROR_DOWN=ERROR;
end;
if (exist('transformfcn'))
    x=PLOT;
    TPLOT=eval(transformfcn);
    x=PLOT+ERROR;
    ERROR_UP=eval(transformfcn)-TPLOT;
    x=PLOT-ERROR;
    ERROR_DOWN=-eval(transformfcn)+TPLOT;
    PLOT=TPLOT;
end;

[numxcat,numxvars]=size(C);
[numsplitcat,numsplitvars]=size(R);
%numlvars=size(D{1,1},2);
glabels=makexlabels(C,xvar_conv);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now format the x-size depending on the grouping structure
l_from=1;
l_to=[];
if (numxvars==1)
    x_coord=C;
else
    x_coord=1;
    for c=2:numxcat
        for gv=1:numxvars
            if(C(c,gv)~=C(c-1,gv))
                x_coord(c,1)=x_coord(c-1)+gap(gv);
                if (gv~=numxvars)
                    l_to(end+1)=c-1;
                    l_from(end+1)=c;
                end;
                break;
            end;
        end;
    end;
    if (strcmp(XCoord,'last'))
        x_coord=C(end,:)';
    end;
end;
l_to(end+1)=numxcat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate scale
xmin = min(x_coord);
xmax = max(x_coord);
dx = (xmax-xmin)/20;if (dx==0) dx=xmax/10;end;
xlims = [(xmin-dx) (xmax+dx)];
if (exist('ERROR','var'))
    ymin = min(min(PLOT-ERROR));
    ymax = max(max(PLOT+ERROR));
else
    ymin = min(min(PLOT));
    ymax = max(max(PLOT));
end;
dy = (ymax-ymin)/20;if (dy==0) dy=ymax/10;end;
ylims = [(ymin-dy) (ymax+dy)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale axis for vertical or horizontal boxes.
% deal with hold on option.
holding=get(gca,'NextPlot');
if (strcmp(holding,'add'))
    xlim_old=get(gca,'XLim');
    ylim_old=get(gca,'YLim');
    xlims=[min(xlims(1),xlim_old(1)) max(xlims(2),xlim_old(2))];
    ylims=[min(ylims(1),ylim_old(1)) max(ylims(2),ylim_old(2))];
else
    cla
    set(gca,'NextPlot','add');
end;
set(gca,'Box','off');
if ~F.flip
    axis([xlims ylims]);
    set(gca,'XTick',x_coord);
    set(gca,'YLabel',text(0,0,'Values'));
    %if (isempty(g)), set(gca,'XLabel',text(0,0,'Column Number')); end
else
    axis([ylims xlims]);
    set(gca,'YTick',lb);
    set(gca,'XLabel',text(0,0,'Values'));
    %if (isempty(g)), set(gca,'YLabel',text(0,0,'Column Number')); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform plotting
for c=1:numsplitcat
    if (~isempty(CAT))
        fields=fieldnames(CAT);
        for f=1:length(fields)
            fiel=getfield(CAT,fields{f});
            size_fiel=length(fiel);
            %F.(fields{f})=fiel{spcat};
            if c>size_fiel
                warning('More categories in split variable than in format structure. Recycling formats!');
            end;
            F=setfield(F,fields{f},fiel{mod(c-1,size_fiel)+1});
        end;
    end;
    for seg=1:length(l_from)
        x=x_coord(l_from(seg):l_to(seg));
        y=PLOT(c,l_from(seg):l_to(seg));
        if (~isempty(errorfcn))
            EU=ERROR_UP(c,l_from(seg):l_to(seg));
            ED=ERROR_DOWN(c,l_from(seg):l_to(seg));
            if (strcmp(F.errorbars,'plusminus'))
                errorbars(x,y',[EU' ED'],'linecolor',F.errorcolor,'linewidth',F.errorwidth,'error_dir','both');
            elseif (strcmp(F.errorbars,'shade'))
                Y=[y+EU fliplr(y-ED)];
                X=[x fliplr(x)];
                h=patch(squeeze(X), squeeze(Y), F.linecolor);
                set (h, 'FaceColor',F.shadecolor,'EdgeColor','none');
            end;
        end;
        h(c)=plot(x,y);
        set(h(c),'Color',F.linecolor,'LineWidth',F.linewidth,'LineStyle',F.linestyle,'Marker',F.markertype,'MarkerSize',F.markersize,'MarkerEdgeColor',F.markercolor,'MarkerFaceColor',F.markerfill);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do legend
if (~isempty(split))
    plotlegend(h,leg,R,split_conv,leglocation);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X-axis labels

%Turn off tick labels and axis label
if (numxvars>1 & XTickLabel==1)
    set(gca, 'XTickLabel','','UserData',numxvars);
    xlabel('');
    ylim = get(gca, 'YLim');

    % Place multi-line text approximately where tick labels belong
    for j=1:numxcat
        ht = text(x_coord(j),ylim(1),glabels{j,1},'HorizontalAlignment','center',...
            'VerticalAlignment','top', 'UserData','xtick');
    end

    % Resize function will position text more accurately
    set(gcf, 'ResizeFcn', sprintf('lineplot(%d)', gcf), ...
        'Interruptible','off', 'PaperPositionMode','auto');
    resizefcn(gcf);
    %set(gca, 'XTickLabel',glabels);


    % Store information for gname: this might be cool!
    % set(gca, 'UserData', {'boxplot' xvisible gorig vert});
elseif (numxvars>1 & XTickLabel==0)
    glabels=makexlabels(C(:,end));
    set(gca,'XTickLabel',glabels);
end;
set(gca,'NextPlot',holding);

function resizefcn(f)
% Adjust figure layout to make sure labels remain visible
h = findobj(f, 'UserData','xtick');
if (isempty(h))
    set(f, 'ResizeFcn', '');
    return;
end
ax = get(f, 'CurrentAxes');
nlines = get(ax, 'UserData');

% Position the axes so that the fake X tick labels have room to display
set(ax, 'Units', 'characters');
p = get(ax, 'Position');
ptop = p(2) + p(4);
if (p(4) < nlines+1.5)
    p(2) = ptop/2;
else
    p(2) = nlines + 1;
end
p(4) = ptop - p(2);
set(ax, 'Position', p);
set(ax, 'Units', 'normalized');

% Position the labels at the proper place
xl = get(gca, 'XLabel');
set(xl, 'Units', 'data');
p = get(xl, 'Position');
ylim = get(gca, 'YLim');
p2 = (p(2)+ylim(1))/2;
for j=1:length(h)
    p = get(h(j), 'Position') ;
    p(2) = p2;
    set(h(j), 'Position', p);
end
