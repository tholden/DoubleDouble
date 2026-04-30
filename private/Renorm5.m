function [ c0, c1, c2, c3 ] = Renorm5( c0, c1, c2, c3, c4 )
    [ c4, c5 ] = Normalize( c3, c4 );
    [ c3, c4 ] = Normalize( c2, c4 );
    [ c2, c3 ] = Normalize( c1, c3 );
    [ c1, c2 ] = Normalize( c0, c2 );
    c0 = c1;

    s0 = c0;
    s1 = c2;

    s2 = zeros( size( s0 ) );
    s3 = zeros( size( s0 ) );

    Select = s1 ~= 0;
    if any( Select, 'all' )
        [ s1( Select ), s2( Select ) ] = Normalize( s1( Select ), c3( Select ) );
        Select2 = Select & ( s2 ~= 0 );
        if any( Select2, 'all' )
            [ s2( Select2 ), s3( Select2 ) ] = Normalize( s2( Select2 ), c4( Select2 ) );
            Select2A = Select2 & ( s3 ~= 0 );
            if any( Select2A, 'all' )
                s3( Select2A ) = s3( Select2A ) + c5( Select2A );
            end
            Select2B = Select2 & ~( s3 ~= 0 );
            if any( Select2B, 'all' )
                s2( Select2B ) = s2( Select2B ) + c5( Select2B );
            end
        end
        Select3 = Select & ~( s2 ~= 0 );
        if any( Select3, 'all' )
            [ s1( Select3 ), s2( Select3 ) ] = Normalize( s1( Select3 ), c4( Select3 ) );
            Select3A = Select3 & ( s2 ~= 0 );
            if any( Select3A, 'all' )
                [ s2( Select3A ), s3( Select3A ) ] = Normalize( s2( Select3A ), c5( Select3A ) );
            end
            Select3B = Select3 & ~( s2 ~= 0 );
            if any( Select3B, 'all' )
                [ s1( Select3B ), s2( Select3B ) ] = Normalize( s1( Select3B ), c5( Select3B ) );
            end
        end
    end

    Select4 = ~Select;
    if any( Select4, 'all' )
        [ s0( Select4 ), s1( Select4 ) ] = Normalize( s0( Select4 ), c3( Select4 ) );
        Select5 = Select4 & ( s1 ~= 0 );
        if any( Select5, 'all' )
            [ s1( Select5 ), s2( Select5 ) ] = Normalize( s1( Select5 ), c4( Select5 ) );
            Select5A = Select5 & ( s2 ~= 0 );
            if any( Select5A, 'all' )
                [ s2( Select5A ), s3( Select5A ) ] = Normalize( s2( Select5A ), c5( Select5A ) );
            end
            Select5B = Select5 & ~( s2 ~= 0 );
            if any( Select5B, 'all' )
                [ s1( Select5B ), s2( Select5B ) ] = Normalize( s1( Select5B ), c5( Select5B ) );
            end
        end
        Select6 = Select4 & ~( s1 ~= 0 );
        if any( Select6, 'all' )
            [ s0( Select6 ), s1( Select6 ) ] = Normalize( s0( Select6 ), c4( Select6 ) );
            Select6A = Select6 & ( s1 ~= 0 );
            if any( Select6A, 'all' )
                [ s1( Select6A ), s2( Select6A ) ] = Normalize( s1( Select6A ), c5( Select6A ) );
            end
            Select6B = Select6 & ~( s1 ~= 0 );
            if any( Select6B, 'all' )
                [ s0( Select6B ), s1( Select6B ) ] = Normalize( s0( Select6B ), c5( Select6B ) );
            end
        end
    end

    c0 = s0; c1 = s1; c2 = s2; c3 = s3;
end
