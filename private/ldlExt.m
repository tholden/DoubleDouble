function [ L, D ] = ldlExt( A, type )
    % Type-agnostic LDL decomposition

    [ m, n ] = size( A );
    assert( m == n, 'Matrix must be square.' );
    L = eye( n, 'like', A );
    D = zeros( 1, n, 'like', A );

    D( 1 ) = A( 1, 1 );
    if n > 1
        L( 2 : n, 1 ) = A( 2 : n, 1 ) ./ D( 1 );
    end

    for j = 2 : n
        idxs = 1 : ( j - 1 );
        t = sum( L( j, idxs ) .* D( idxs ) .* conj( L( j, idxs ) ) );
        D( j ) = A( j, j ) - t;

        if j < n
            jdxs = ( j + 1 ) : n;
            tt = sum( L( jdxs, idxs ) .* ( D( idxs ) .* conj( L( j, idxs ) ) ), 2 );
            L( jdxs, j ) = ( A( jdxs, j ) - tt ) ./ D( j );
        end
    end

    if nargin < 2 || ~strcmp( type, 'vector_d' )
        D = diag( D );
    else
        D = D.';
    end

end
