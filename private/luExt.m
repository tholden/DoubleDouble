function [ L, U, P ] = luExt( A )
    % Type-agnostic LU decomposition with partial pivoting
    [ m, n ] = size( A );
    P = eye( m );
    L = eye( m, n, 'like', A );
    U = A;
    for k = 1 : min( m, n )
        [ ~, max_idx ] = max( abs( U( k : m, k ) ) );
        max_idx = max_idx + k - 1;
        
        if max_idx ~= k
            U( [ k, max_idx ], : ) = U( [ max_idx, k ], : );
            P( [ k, max_idx ], : ) = P( [ max_idx, k ], : );
            if k > 1
                L( [ k, max_idx ], 1 : k - 1 ) = L( [ max_idx, k ], 1 : k - 1 );
            end
        end
        
        if U( k, k ) ~= 0
            idx = k + 1 : m;
            L( idx, k ) = U( idx, k ) ./ U( k, k );
            U( idx, k : n ) = U( idx, k : n ) - L( idx, k ) * U( k, k : n );
        end
    end
end
