function [ Q, R ] = qrExt( A )
    % Type-agnostic QR decomposition (Modified Gram-Schmidt)
    [ m, n ] = size( A );
    Q = zeros( m, n, 'like', A );
    R = zeros( n, n, 'like', A );
    
    for j = 1 : n
        v = A( :, j );
        for i = 1 : j - 1
            R( i, j ) = sum( conj( Q( :, i ) ) .* v );
            v = v - R( i, j ) * Q( :, i );
        end
        R( j, j ) = norm( v );
        if R( j, j ) > 0
            Q( :, j ) = v ./ R( j, j );
        else
            Q( :, j ) = v;
        end
    end
end
