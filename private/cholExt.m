function [ v, p ] = cholExt( v, type )
    [ v, d ] = ldlExt( v, 'vector_d' );
    v = v .* sqrt( d.' );
    if any( d < 0 )
        p = 1;
    else
        p = 0;
    end
    if nargin < 2 || strcmp( type, 'upper' )
        v = v.';
    end
end
