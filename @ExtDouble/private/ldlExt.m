function [ L, D ] = ldlExt( A, type ) % Type-agnostic LDL decomposition

    [ m, n ] = size( A );
    assert( m == n, 'Matrix must be square.' );
    L = eye( n, 'like', A );
    D = zeros( 1, n, 'like', A );

    D = D.Assign( A.Index( 1, 1 ), 1 );
    if n > 1
        L = L.Assign( A.Index( 2 : n, 1 ) ./ D.Index( 1 ), 2 : n, 1 );
    end

    for j = 2 : n
        idxs = 1 : ( j - 1 );
        t = sum( L.Index( j, idxs ) .* D.Index( idxs ) .* conj( L.Index( j, idxs ) ) );
        D = D.Assign( A.Index( j, j ) - t, j );

        if j < n
            jdxs = ( j + 1 ) : n;
            tt = sum( L.Index( jdxs, idxs ) .* ( D.Index( idxs ) .* conj( L.Index( j, idxs ) ) ), 2 );
            L = L.Assign( ( A.Index( jdxs, j ) - tt ) ./ D.Index( j ), jdxs, j );
        end
    end

    if nargin < 2 || ~strcmp( type, 'vector_d' )
        D = diag( D );
    else
        D = D.';
    end

end
