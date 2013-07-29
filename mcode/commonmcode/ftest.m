function [h, p] = ftest(x, y, alpha)
    xbar = sum(x) / length(x);
    ybar = sum(y) / length(y);
    sq1 = sum((x - xbar) .^ 2) / (length(x) - 1);
    sq2 = sum((y - ybar) .^ 2) / (length(y) - 1);
    F = sq1 / sq2;
    Fc1 = finv(alpha / 2, length(x) - 1, length(y) - 1);
    Fc2 = finv(1 - alpha / 2, length(x) -1, length(y) - 1);
    % Two-tailed test
    if (F > Fc2 | F < Fc1);
        h = 1;  % null hypothesis rejected
    else
        h = 0;
    end
    % calculate the p-value
    p = fcdf(F, length(x) - 1, length(y) - 1);    
return