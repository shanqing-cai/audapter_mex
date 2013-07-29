function output = mva(input, width, varargin)    % width must be equal to or greater than 3
    v = input;
    if(size(v, 1) > size(v, 2))
        v = v';
    end
    if (isempty(input))
        output = input;
        return
    end
    
    % Padding
    v = [ones(1, width - 1) * v(1), v, ones(1, width - 1) * v(end)];
    
    kern = [1:1:ceil(width / 2), ceil(width/2) - 1:-1:1];   % A triangular window whose width is always odd number
    kern = kern / sum(kern);
    if (~isempty(findStringInCell(varargin,'hamming')) | ~isempty(findStringInCell(varargin,'Hamming')))
        kern=hamming(width);
        if (size(kern,1)>size(kern,2))
            kern=kern';
        end
        kern=kern/sum(kern);
    end    
    output = conv(v, kern);
    
    while(length(output) > length(input))
        output = output(2:end - 1);     % The output is ensured to have the same length as the input
    end
        
return