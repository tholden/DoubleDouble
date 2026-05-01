function [ c0, c1, c2, c3 ] = Renorm4( c0, c1, c2, c3 )
    [ c3, c4 ] = EDNormalize( c2, c3 );
    [ c2, c3 ] = EDNormalize( c1, c3 );
    [ c1, c2 ] = EDNormalize( c0, c2 );
    c0 = c1;

    s0 = c0;
    s1 = c2;

    s2 = zeros( size( s0 ) );
    s3 = zeros( size( s0 ) );

    Select = s1 ~= 0;
    if any( Select, 'all' )
        [ s1( Select ), s2( Select ) ] = EDNormalize( s1( Select ), c3( Select ) );
        Select2 = Select & ( s2 ~= 0 );
        if any( Select2, 'all' )
            [ s2( Select2 ), s3( Select2 ) ] = EDNormalize( s2( Select2 ), c4( Select2 ) );
        end
        Select3 = Select & ~( s2 ~= 0 );
        if any( Select3, 'all' )
            [ s1( Select3 ), s2( Select3 ) ] = EDNormalize( s1( Select3 ), c4( Select3 ) );
        end
    end

    Select4 = ~Select;
    if any( Select4, 'all' )
        [ s0( Select4 ), s1( Select4 ) ] = EDNormalize( s0( Select4 ), c3( Select4 ) );
        Select5 = Select4 & ( s1 ~= 0 );
        if any( Select5, 'all' )
            [ s1( Select5 ), s2( Select5 ) ] = EDNormalize( s1( Select5 ), c4( Select5 ) );
        end
        Select6 = Select4 & ~( s1 ~= 0 );
        if any( Select6, 'all' )
            [ s0( Select6 ), s1( Select6 ) ] = EDNormalize( s0( Select6 ), c4( Select6 ) );
        end
    end

    c0 = s0; c1 = s1; c2 = s2; c3 = s3;
end
