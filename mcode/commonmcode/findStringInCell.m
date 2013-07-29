function index =  findStringInCell(theCell, theString, varargin)
    index = [];
    
    inverse=0;
    if (nargin>2)
    	if (isequal(varargin{1},'inv') | isequal(varargin{1},'inverse'))
            inverse=1;
        end
    end
    for n = 1:length(theCell)
        if (~isstruct(theCell{n}))
            if (~inverse)
                if (strcmp(theCell{n}, theString))
                    index = [index, n];
                end
            else
                if (~strcmp(theCell{n}, theString))
                    index = [index, n];
                end
            end
        end
    end
return