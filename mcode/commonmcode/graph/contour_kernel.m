function [X,Y,Z]=contour_kernel(x,y,z,varargin) 
% Makes a countour plot for a kernel estimate 
% either of the density of x,y pairs, 
% or of the mean of the z-variable
% 
% function contour_pivot(x,y,z,stats) 
% x: column vector of X values
% y: column vector of Y values
% z: column vector of Z values: if empty it plots density estimate
% Varargin
%   'stats',{'mean','count'}: default mean
%   'numcat',[x,y]      : number of categories in X and Y directions(resolution)
%   'vertlines',a       : Indices for vertical lines on top of countour plot 
%   'sigma',[SDx,SDy]   : SD of smoothing kernel in x and y directions  
%   'subset',logical    : Subset index 
if (nargin<3 | isempty (z))
    stats='count';
else
    stats='mean';
end;
numcat=[];
sigma=[];
vertlines=[];
subset=(x==x);
% Deal with the varargin's 
c=1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'stats','numcat','vertlines','sigma','subset'}
            eval([varargin{c} '=varargin{c+1};']);
            c=c+2;
    end;        
end;
if (nargin<3)
    z=[];
end;

x=x(find(subset));
y=y(find(subset));
if (~isempty(z))
    z=z(find(subset));
end;
if (isempty(numcat))
    numcat=[floor(length(x)/10) floor(length(x)/10)];
end;
if (isempty(sigma))
    sigma=[(max(x)-min(x))/numcat(1) (max(y)-min(y))/numcat(2)]*2.5;
end;
stepx=(max(x)-min(x))/numcat(1);
stepy=(max(y)-min(y))/numcat(2);
[X,Y]=meshgrid([min(x)+stepx/2:stepx:max(x)-stepx/2],[min(y)+stepy/2:stepy:max(y)-stepy/2]);
if (strcmp(stats,'mean'))
    if (all(sigma==0)) 
        [N,Z]=cat2d(x,y,z,X,Y);
    else 
        [N,Z]=kernel_est2d(x,y,z,X,Y,sigma);
    end;
    [c,h]=contourf(X,Y,Z);
else
    N=kernel_est2d(x,y,z,X,Y,sigma);
    [c,h]=contourf(X,Y,N);
end;
clabel(c,h);
if (~isempty(vertlines))
    drawlines(find(vertlines(find(subset))),'k');
end;
