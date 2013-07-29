function [F,R,C]=contour_pivot(X,Y,Z,stats,varargin) 
% makes a countour plot of some statistics from a data-columns
% function contour_pivot(Y,X,Z,stats) 
% X: column vector of X values
% Y: column vector of Y values
% Z: column vector of Z values
% stats: name of the function applied to each collection of Z-values for a particular X,Y
c=1;PV={};
subset=ones(length(stats),1);
contours=[];
vararginoptions(varargin,{'subset','contours'});

[F,R,C]=pivottable(Y,X,Z,stats,'subset',subset);
[x,y]=meshgrid(C,R);
if(strcmp(stats,'length'))
    F(isnan(F))=0;
end;
if (isempty(contours)) 
    [c,h]=contourf(x,y,F);
else 
    [c,h]=contourf(x,y,F,contours);
end;    
% colorbar;
clabel(c,h);