function v = IsEqualWithExpansion( a, b, varargin )
    v = a == b;
    v = all( v(:) );
    if nargin > 2
        for i = 1 : length( varargin )
            if ~v
                break
            end
            v = v && IsEqualWithExpansion( a, varargin{ i } );
        end
    end
end
