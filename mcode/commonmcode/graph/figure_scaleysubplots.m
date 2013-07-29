function figure_scaleysubplots(subplt)
% figure_scaleysubplots(subplt)
% Make the y-scale on all subplots subplt the same 
% Goes after minimal min and maximal max
% if subplt not given, does it for all
h=get(gcf,'Children');
if nargin==0
    subplt=1:length(h)
end;
for i=subplt
    y(i,:)=get(h(i),'YLim');
end;
ylim=[min(y(:,1)) max(y(:,2))];
for i=subplt
    set(h(i),'YLim',ylim);
end;

    
