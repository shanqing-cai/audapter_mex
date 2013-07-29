function is = isfile(fileName)
    d = dir(fileName);
    if (length(d) == 0)
        is = 0;
    else
        is = 1;
    end
return