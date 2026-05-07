function [ x0, x1, x2, x3 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % IEEE-accurate QD+QD addition ( cf. qd_real::ieee_add in QD library ).
    % Sorts all 8 components by decreasing magnitude, then cascades via
    % the sloppy_add accumulation structure for vectorized operation.

    if ~isreal( a0 ) || ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( a3 ) || ...
            ~isreal( b0 ) || ~isreal( b1 ) || ~isreal( b2 ) || ~isreal( b3 )
        [ x0r, x1r, x2r, x3r ] = QDPlusQD( real( a0 ), real( a1 ), real( a2 ), real( a3 ),  ...
            real( b0 ), real( b1 ), real( b2 ), real( b3 ) );
        [ x0i, x1i, x2i, x3i ] = QDPlusQD( imag( a0 ), imag( a1 ), imag( a2 ), imag( a3 ),  ...
            imag( b0 ), imag( b1 ), imag( b2 ), imag( b3 ) );
        x0 = complex( x0r, x0i );
        x1 = complex( x1r, x1i );
        x2 = complex( x2r, x2i );
        x3 = complex( x3r, x3i );
        return
    end

    [ a0, a1, a2, a3, b0, b1, b2, b3 ] = ED.ExpandSingleton( a0, a1, a2, a3, b0, b1, b2, b3 );
    CatDim = ndims( a0 ) + 1;
    z = cat( CatDim, a0, b0, a1, b1, a2, b2, a3, b3 );
    z = sort( z, CatDim, 'descend', 'ComparisonMethod', 'abs' );

    sz = size( z );
    zFlat = reshape( z, [], sz( end ) );

    u = zFlat( :, 1 );
    v = zFlat( :, 2 );

    % quick_two_sum(u, v, v)
    s_uv = u + v;
    v = v - (s_uv - u);
    u = s_uv;

    x1 = zeros( size( u ), 'like', u );
    x2 = zeros( size( u ), 'like', u );
    x3 = zeros( size( u ), 'like', u );
    x4 = zeros( size( u ), 'like', u );

    k = ones( size( u ), 'like', u );

    for i = 3 : 8
        t = zFlat( :, i );
        Active = k <= 4;

        if any( Active )
            va = v( Active );
            ta = t( Active );
            ua = u( Active );

            % s = two_sum(b, c, b);
            s1 = va + ta;
            bb1 = s1 - va;
            v_new = (va - (s1 - bb1)) + (ta - bb1);

            % s = two_sum(a, s, a);
            s2 = ua + s1;
            bb2 = s2 - ua;
            u_new = (ua - (s2 - bb2)) + (s1 - bb2);

            za = u_new ~= 0;
            zb = v_new ~= 0;
            emit_act = za & zb;

            u_final = s2;
            v_final = v_new;

            cond_not_zb = ~zb;
            v_final( cond_not_zb ) = u_new( cond_not_zb );

            u_final( emit_act ) = u_new( emit_act );
            v_final( emit_act ) = v_new( emit_act );

            u( Active ) = u_final;
            v( Active ) = v_final;

            emit = false( size( u ) );
            emit( Active ) = emit_act;

            s_full = zeros( size( u ), 'like', u );
            s_full( Active ) = s2;

            old_k = k;

            idx1 = emit & (old_k == 1); x1( idx1 ) = s_full( idx1 ); k( idx1 ) = 2;
            idx2 = emit & (old_k == 2); x2( idx2 ) = s_full( idx2 ); k( idx2 ) = 3;
            idx3 = emit & (old_k == 3); x3( idx3 ) = s_full( idx3 ); k( idx3 ) = 4;
            idx4 = emit & (old_k == 4); x4( idx4 ) = s_full( idx4 ); k( idx4 ) = 5;
        end

        inActive = ~Active;
        if any( inActive )
            x4( inActive ) = x4( inActive ) + t( inActive );
        end
    end

    idx1 = k == 1; x1( idx1 ) = u( idx1 ); x2( idx1 ) = v( idx1 );
    idx2 = k == 2; x2( idx2 ) = u( idx2 ); x3( idx2 ) = v( idx2 );
    idx3 = k == 3; x3( idx3 ) = u( idx3 ); x4( idx3 ) = v( idx3 );
    idx4 = k == 4; x4( idx4 ) = u( idx4 );

    [ x0, x1, x2, x3 ] = QDNormalize( x1, x2, x3, x4 );

    x0 = reshape( x0, sz( 1 : ( end - 1 ) ) );
    x1 = reshape( x1, sz( 1 : ( end - 1 ) ) );
    x2 = reshape( x2, sz( 1 : ( end - 1 ) ) );
    x3 = reshape( x3, sz( 1 : ( end - 1 ) ) );

end
