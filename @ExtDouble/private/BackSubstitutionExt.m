function v = BackSubstitutionExt( v, U ) % For upper triangular U, x = BackSubstitutionExt( b, U ) solves U*x = b.

    [ m, n ] = size( U );
    mn = min( m, n );
    [ vm, vn ] = size( v );
    if vm < n
        v = [ v; zeros( n - vm, vn, 'like', v ) ];
    elseif vm > n
        v = v.Index( 1 : n, ':' );
    end

    v = v.Assign( v.Index( mn, ':' ) ./ U.Index( mn, mn ), mn, ':' );
    for k = ( mn - 1 ) : -1 : 1
        j = ( k + 1 ) : n;
        t = sum( v.Index( j, ':' ) .* U.Index( k, j ).', 1 );
        v = v.Assign( ( v.Index( k, ':' ) - t ) ./ U.Index( k, k ), k, ':' );
    end

end
