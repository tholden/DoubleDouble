function w = convExt( u, v )

    RowVector = size( u, 1 ) == 1 && size( v, 1 ) == 1;

    u = u.Vec();
    v = v.Vec();

    M = size( u, 1 );
    N = size( v, 1 );

    K = M + N - 1;

    w = zeros( K, 1, 'like', u );

    for k = 1 : K

        j = max( 1, k + 1 - N ) : min( k, M );
        i = k - j + 1;

        wk = dot( u.Index( j ), v.Index( i ) );
        w = w.Assign( wk, k );

    end

    if RowVector
        w = w.';
    end

end
