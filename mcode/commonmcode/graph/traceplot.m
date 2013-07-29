function PLOT=traceplot(t,data,varargin);
% function T=traceplot(time,rawtraces,varargin);
% makes a plot to compare the mean of temporal time series between conditions
% t: is a 1xT vector of times 
% data: is a NxT matrix of time series 
% VARARGIN:
%   'subset',index: restrict the plotting to a subset of the data 
%   'split',var: Plot different lines for each category of var 
%   'leg','auto'/{texts} : add a legend to the figure
%   'errorfcn','std'/'stderr': Add shapded error bar area behind trace
% formating options: 
%   'linecolor',[r g b]: Line color 
%   'linestyle',{'-',':',...}: Line Style option
%   'linewidth',x: line width 
% v 1.1 Displaying one timeseries per category is possible now
%       recycling of format in split categories 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defaults and deal with options 
leg={};
plotfcn='nanmean';
errorfcn='';
subset=[];
split=[];
linewidth=1;
linecolor=[0 0 0];
linestyle='-';
patchcolor=[0.8 0.8 0.8];
transp=0.3;
XLim=[];
YLim=[];
options={'subset','split','leg','plotfcn','errorfcn','linewidth','CAT','linecolor','linestyle','XLim','YLim'};
flags={};
CAT.linecolor={[1 0 0],[0 0 1],[0 1 0],[0 0 0],[1 1 0]};
CAT.patchcolor={[1 0 0],[0 0 1],[0 1 0],[0 0 0],[1 1 0]};
vararginoptions(varargin,options,flags);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check proper sizing of Input arguments 
[a,T]=size(t);
if (a>1) 
    t=t';
    [a,T]=size(t);
end;
if (a>1) error('Time vector must be a vector, not matrix'); end;
[N,T1]=size(data);
if (T~=T) error('Data must be a NxT matrix'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with subset and splitby arguments 
if (isempty(split))
    split=ones(N,1);
end;
[split,split_conv]=fac2int(split);
if (~isempty(subset))
    data=data(subset,:);
    split=split(subset,:);
end;
A=pidata(split,data);
numcats=size(A,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the plotting 
for c=1:numcats
    if (~isempty(CAT))
        fields=fieldnames(CAT);
        for f=1:length(fields)
            fiel=getfield(CAT,fields{f});
            a=mod(c-1,length(fiel))+1;
            if (c>length(fiel))
                warning('reusing format CAT');
            end;
            eval([fields{f} '=fiel{a};']);
        end;
    end;
    if (size(A{c,2},1)==1) 
        PLOT(c,:)=A{c,2};
    else
        PLOT(c,:)=fcneval(plotfcn,(A{c,2}));
    end;
    if (~isempty(errorfcn))
        ERR(c,:)=fcneval(errorfcn,A{c,2});
        p(c)=plotshade(t,PLOT(c,:),ERR(c,:),'patchcolor',patchcolor,'transp',transp);hold on;
    end;
    h(c)=plot(t,PLOT(c,:));hold on;
    set(h(c),'LineWidth',linewidth,'Color',linecolor,'LineStyle',linestyle)
end;
hold off;
set(gca,'Box','off');
if (~isempty(XLim))
    set(gca,'XLim',XLim);
end;
if (~isempty(YLim))
    set(gca,'YLim',YLim);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do legend 
if (~isempty(split))
    R=vertcat(A{:,1});
    plotlegend(h,leg,R,split_conv);
end;
