function [ c0, c1, c2, c3 ] = QDNormalize( c0, c1, c2, c3 )

    if ~isreal( c0 ) || ~isreal( c1 ) || ~isreal( c2 ) || ~isreal( c3 )
        [ c0r, c1r, c2r, c3r ] = QDNormalize( real( c0 ), real( c1 ), real( c2 ), real( c3 ) );
        [ c0i, c1i, c2i, c3i ] = QDNormalize( imag( c0 ), imag( c1 ), imag( c2 ), imag( c3 ) );
        c0 = complex( c0r, c0i );
        c1 = complex( c1r, c1i );
        c2 = complex( c2r, c2i );
        c3 = complex( c3r, c3i );
        return
    end

    c0 = Normalize( c0 );
    c1 = Normalize( c1 );
    c2 = Normalize( c2 );
    c3 = Normalize( c3 );

    [ c3, c4 ] = DDNormalize( c2, c3 );
    [ c2, c3 ] = DDNormalize( c1, c3 );
    [ c1, c2 ] = DDNormalize( c0, c2 );

    c0 = c1;

    s0 = c0;
    s1 = c2;

    s2 = zeros( size( s0 ), 'like', s0 );
    s3 = zeros( size( s0 ), 'like', s0 );

    Select = s1 ~= 0;
    if any( Select, 'all' )
        [ s1( Select ), s2( Select ) ] = DDNormalize( s1( Select ), c3( Select ) );
        Select2 = Select & ( s2 ~= 0 );
        if any( Select2, 'all' )
            [ s2( Select2 ), s3( Select2 ) ] = DDNormalize( s2( Select2 ), c4( Select2 ) );
        end
        Select3 = Select & ~( s2 ~= 0 );
        if any( Select3, 'all' )
            [ s1( Select3 ), s2( Select3 ) ] = DDNormalize( s1( Select3 ), c4( Select3 ) );
        end
    end

    Select4 = ~Select;
    if any( Select4, 'all' )
        [ s0( Select4 ), s1( Select4 ) ] = DDNormalize( s0( Select4 ), c3( Select4 ) );
        Select5 = Select4 & ( s1 ~= 0 );
        if any( Select5, 'all' )
            [ s1( Select5 ), s2( Select5 ) ] = DDNormalize( s1( Select5 ), c4( Select5 ) );
        end
        Select6 = Select4 & ~( s1 ~= 0 );
        if any( Select6, 'all' )
            [ s0( Select6 ), s1( Select6 ) ] = DDNormalize( s0( Select6 ), c4( Select6 ) );
        end
    end

    c0 = s0;
    c1 = s1;
    c2 = s2;
    c3 = s3;

end
