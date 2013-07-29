function figure_scalesubplots()
% function figure_scalesubplots()
%   scales all subplots in a figure to the same axis 
% see also figure_scalexsubplots, figure_scaleysubplots
% Joern Diedrichsen, jdiedric@jhu.edu
% v.1.0 9/13/05
h=get(gcf,'Children');
for i=1:length(h);
    x(i,:)=get(h(i),'XLim');
    y(i,:)=get(h(i),'YLim');
end;
xlim=[min(x(:,1)) max(x(:,2))];
ylim=[min(y(:,1)) max(y(:,2))];
for i=1:length(h);
    set(h(i),'XLim',xlim);
    set(h(i),'YLim',ylim);
end;

    
