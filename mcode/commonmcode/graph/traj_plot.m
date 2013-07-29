function traj_plot(X,Y,varargin);
%  function traj_plot(x,y,varargin)
% Description
%  varargin:
%   Format options (for all symbols) 
%       'linecolor': 
%       'markercolor':
%       'markersize':
%       'markersize': 
%   Other commands:
%       'CAT': Stucture with fields with one entry per category (split by) (all fields are 
%                 cell arrays with length cat
%       'subset': plots only a subset of the data
% Joern Diedrichsen jdiedric@jhu.edu
[tracelength,numtraces]=size(X);
[tracelengthy,numtracesy]=size(Y);
if (tracelength~=tracelengthy | numtraces~=numtracesy)
    error('X and Y have to have the same size\n');
end;
% Set defaults for all plots 
F.markertype =NaN;
F.markercolor=[0 0 0];
F.markerfill=[1 1 1];
F.markersize=4;
F.linecolor=[0 0 0];
F.linestyle='-';
F.linewidth=1;
draworig=0;
style='single';
CAT.markercolor={[1 0 0],[0 0 0],[1 0 0],[0 0 0]};
CAT.linecolor={[0 0 0],[1 0 0],[0 1 0],[0 0 1]};

% Deal with the varargin's 
c=1;
split=[];
while(c<=length(varargin))
    switch(varargin{c})
        case 'split'
            split=varargin{c+1};c=c+2;
        case 'markersize'
            F.markersize=varargin{c+1};c=c+2;
        case 'markertype'
            F.markertype=varargin{c+1};c=c+2;
        case 'linecolor'
            F.linecolor=varargin{c+1};c=c+2;
        case 'CAT'
            CAT=varargin{c+1};c=c+2;
        case 'subset'
            goodindx=find(varargin{c+1});
            c=c+2;
            Y=Y(:,goodindx);
            X=X(:,goodindx);
            if (~isempty(split))            
                split=split(goodindx,:);
            end;
        otherwise
            error('Unknown option\n');
    end;
end;

% Calculate scale 
xmin = min(min(X));
xmax = max(max(X));
dx = (xmax-xmin)/20;
xlims = [(xmin-dx) (xmax+dx)];

ymin = min(min(Y));
ymax = max(max(Y));
dy = (ymax-ymin)/20;
ylims = [(ymin-dy) (ymax+dy)];


% deal with split 
set(gca,'NextPlot','add','Box','off');
if (~isempty(split))
    [numobserv,numsplitvar] 
    D=pidata(split,[1:num_traces]');
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
        do_traj_plot(X(:,D{c,2}),Y(:,D{c,2}),F,style);
    end;
else
    do_traj_plot(X,Y,F,style);
end;
axis([xlims ylims]);
if (draworig)
    draworigin;
end;
set(gca,'NextPlot','replace');

function do_traj_plot(X,Y,F,style);
[dummy,numtraces]=size(X);
if strcmp(style,'single')
    for n=1:numtraces
        if (~isnan(F.markertype))
            set(plot(X(:,n),Y(:,n)),'Marker',F.markertype,'MarkerSize',F.markersize,...
                'MarkerEdgeColor',F.markercolor,'MarkerFaceColor',F.markerfill);
        end;
        if (~isnan(F.linecolor))
            set(plot(X(:,n),Y(:,n),'k'),'Color',F.linecolor,'LineStyle',F.linestyle,'LineWidth',F.linewidth);
        end;        
        
    end;
end;

