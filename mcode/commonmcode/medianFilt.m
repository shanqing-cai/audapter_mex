function y=medianFilt(x,k)
y=zeros(size(x));

h=floor(k/2);
for n=1:length(x)
    y(n)=median(x(max(n-h,1):min(n+h,length(x))));
end

return