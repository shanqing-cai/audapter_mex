function cv=coeffOfVar(x)
idxnn=find(~isnan(x));
x=x(idxnn);
cv=std(x)/mean(x);
return