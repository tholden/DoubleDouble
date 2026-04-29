function v = ForwardEliminationExt( v, L )
    % For lower triangular L, x = ForwardEliminationExt( b, L ) solves L*x = b.

    [ m, n ] = size( L );
    mn = min( m, n );
    [ vm, vn ] = size( v );
    if vm < n
        v = [ v; zeros( n - vm, vn, 'like', v ) ];
    elseif vm > n
        v = v( 1 : n, : );
    end

    v( 1, : ) = v( 1, : ) ./ L( 1, 1 );
    for k = 2 : mn
        j = 1 : ( k - 1 );
        t = sum( v( j, : ) .* L( k, j ).', 1 );
        v( k, : ) = ( v( k, : ) - t ) ./ L( k, k );
    end

end
