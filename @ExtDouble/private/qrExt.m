function [ Q, R ] = qrExt( A ) % Type-agnostic QR decomposition (Modified Gram-Schmidt)

    [ m, n ] = size( A );
    Q = zeros( m, n, 'like', A );
    R = zeros( n, n, 'like', A );

    for j = 1 : n
        v = A.Index( ':', j );
        for i = 1 : ( j - 1 )
            R = R.Assign( sum( conj( Q.Index( ':', i ) ) .* v ), i, j );
            v = v - R.Index( i, j ) * Q.Index( ':', i );
        end
        R = R.Assign( norm( v ), j, j );
        if R.Index( j, j ) > 0
            Q = Q.Assign( v ./ R.Index( j, j ), ':', j );
        else
            Q = Q.Assign( v, ':', j );
        end
    end

end
