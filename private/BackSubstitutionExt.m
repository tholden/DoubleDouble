function v = BackSubstitutionExt( v, U )
    % For upper triangular U, x = BackSubstitutionExt( b, U ) solves U*x = b.

    [ m, n ] = size( U );
    mn = min( m, n );
    [ vm, vn ] = size( v );
    if vm < n
        v = [ v; zeros( n - vm, vn, 'like', v ) ];
    elseif vm > n
        v = v( 1 : n, : );
    end

    v( mn, : ) = v( mn, : ) ./ U( mn, mn );
    for k = ( mn - 1 ) : -1 : 1
        j = ( k + 1 ) : n;
        t = sum( v( j, : ) .* U( k, j ).', 1 );
        v( k, : ) = ( v( k, : ) - t ) ./ U( k, k );
    end

end
