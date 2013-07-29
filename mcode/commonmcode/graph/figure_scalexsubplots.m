function figure_scalexsubplots(subplt)
h=get(gcf,'Children');
if nargin==0
    subplt=1:length(h)
end;
for i=subplt
    x(i,:)=get(h(i),'XLim');
end;
xlim=[min(x(:,1)) max(x(:,2))];
for i=subplt
    set(h(i),'XLim',xlim);
end;

    
