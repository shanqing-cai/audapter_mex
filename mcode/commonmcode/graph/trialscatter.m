function trial_scatter(tn,y,varargin);
% function trial_scatter(tn,y,varargin);
% tn: vector of trial numbers 
% y: plotted measure
%   'split'
% Variable subgroups can be plotted 
% colors: {'k.','r.','ko','ro'}
format={'k.','r.','ko','ro'};
t=[1:length(tn)]';
first=find(tn==1);
option={};BL=[];
c=1;
CAT.markertype= {'.','.','.','.','.','.','.','.'};
CAT.markercolor={[0 0 0],[1 0 0],[0 0 1],[0 1 0],[1 1 0],[0 1 1],[1 0 1],[0.8 0.8 0.8]};

while c<=length(varargin)
    if (strcmp(varargin{c},'blocklabel'))
        BL=varargin{c+1};c=c+2;
    else
        option{end+1}=varargin{c};c=c+1;
    end;    
end;
scatterplot(t,y,'markertype','.','markersize',10,'CAT',CAT,option{:});
drawline(first);
if (~isempty(BL));
    YLIM=get(gca,'YLim');
    y=(YLIM(2)-YLIM(1))*0.85+YLIM(1);
    for i=1:length(first)
        text(first(i),y,num2str(BL(first(i)),'%2.2d'));
    end;
end;