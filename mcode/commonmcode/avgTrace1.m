function r=avgTrace1(xs)
    len=0;
    for n=1:length(xs)
        if (length(xs{n})>len)
            len=length(xs{n});
        end
    end
    
    r=nan(len,3);   % avg(x), std(x), n(x);
    bins=cell(1,len);
    for n=1:len
        for m=1:length(xs)
            if (length(xs{m})>=n)
                bins{n}=[bins{n},xs{m}(n)];
            end
        end
    end
    for n=1:len
        idxnn=find(~isnan(bins{n}));
        r(n,1)=mean(bins{n}(idxnn));
        r(n,2)=std(bins{n}(idxnn));
        r(n,3)=length(bins{n}(idxnn));
    end
return