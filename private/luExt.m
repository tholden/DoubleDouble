function [ v, U, p ] = luExt( v, type )
    [ m, n ] = size( v );
    p = 1 : m;

    for k = 1 : min( m, n )

        % Find index of largest element below diagonal in k-th column
        [ ~, midx ] = max( abs( v( k:m, k ) ) );
        midx = midx + k - 1;

        % Skip elimination if column is zero
        if v( midx, k ) ~= 0

            % Swap pivot row
            if midx ~= k
                v( [ k midx ], : ) = v( [ midx k ], : );
                p( [ k midx ] ) = p( [ midx k ] );
            end

            % Compute multipliers
            i = ( k + 1 ) : m;
            v( i, k ) = v( i, k ) ./ v( k, k );

            % Update the remainder of the matrix
            j = ( k + 1 ) : n;
            v( i, j ) = v( i, j ) - v( i, k ) .* v( k, j );
        end
    end

    if nargout > 1
        % Separate result
        L = tril( v, -1 ) + eye( m, n, 'like', v );
        U = triu( v );
        if n > m
            L = L( :, 1:m );
        elseif n < m
            U = U( 1:n, : );
        end
        v = L;

        if nargout > 2
            if nargin < 2 || ~strcmp( type, 'vector' )
                pp = eye( m );
                pp = pp( p, : );
                p = pp;
            end
        else
            invp( p ) = 1 : m;
            v = v( invp, : );
        end
    end
end
