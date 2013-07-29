function p=chi2test(tab1)
% tab1 is a Mx2 table
    chi2=0;
    M=size(tab1,1);
    N=sum(sum(tab1));
    for n=1:size(tab1,1)
        for m=1:size(tab1,2)
            chi2=chi2+(tab1(n,m)-sum(tab1(n,:))*sum(tab1(:,m))/N)^2 /(sum(tab1(n,:))*sum(tab1(:,m))/N);
        end
    end    
    p=1-chi2cdf(chi2,M);   
end
