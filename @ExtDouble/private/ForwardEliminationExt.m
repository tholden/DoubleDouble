function v = ForwardEliminationExt( v, L ) % For lower triangular L, x = ForwardEliminationExt( b, L ) solves L*x = b.

    [ m, n ] = size( L );
    mn = min( m, n );
    [ vm, vn ] = size( v );
    if vm < n
        v = [ v; zeros( n - vm, vn, 'like', v ) ];
    elseif vm > n
        v = v.Index( 1 : n, ':' );
    end

    v = v.Assign( v.Index( 1, ':' ) ./ L.Index( 1, 1 ), 1, ':' );
    for k = 2 : mn
        j = 1 : ( k - 1 );
        t = sum( v.Index( j, ':' ) .* L.Index( k, j ).', 1 );
        v = v.Assign( ( v.Index( k, ':' ) - t ) ./ L.Index( k, k ), k, ':' );
    end

end
