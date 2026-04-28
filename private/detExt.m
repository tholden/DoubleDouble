function v = detExt( v )
    [ m, n ] = size( v );
    if m ~= n
        throw( MException( 'MATLAB:square', 'Matrix must be square.' ) );
    end
    [ ~, u, P ] = luExt( v );
    DetP = det( P );
    if DetP > 0
        v = prod( diag( u ) );
    elseif DetP < 0
        v = -prod( diag( u ) );
    else
        v = NaN( 'like', v );
    end
end
