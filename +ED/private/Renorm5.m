function [ c0, c1, c2, c3 ] = Renorm5( c0, c1, c2, c3, c4 )
    [ c3, c4 ] = DDNormalize( c3, c4 );
    [ c2, c3 ] = DDNormalize( c2, c3 );
    [ c1, c2 ] = DDNormalize( c1, c2 );
    [ c0, c1 ] = DDNormalize( c0, c1 );

    s0 = c0;
    s1 = c1;

    s2 = zeros( size( s0 ), 'like', s0 );
    s3 = zeros( size( s0 ), 'like', s0 );

    Select = s1 ~= 0;
    if any( Select, 'all' )
        [ s1( Select ), s2( Select ) ] = DDNormalize( s1( Select ), c2( Select ) );
        Select2 = Select & ( s2 ~= 0 );
        if any( Select2, 'all' )
            [ s2( Select2 ), s3( Select2 ) ] = DDNormalize( s2( Select2 ), c3( Select2 ) );
            Select2A = Select2 & ( s3 ~= 0 );
            if any( Select2A, 'all' )
                s3( Select2A ) = s3( Select2A ) + c4( Select2A );
            end
            Select2B = Select2 & ~( s3 ~= 0 );
            if any( Select2B, 'all' )
                s2( Select2B ) = s2( Select2B ) + c4( Select2B );
            end
        end
        Select3 = Select & ~( s2 ~= 0 );
        if any( Select3, 'all' )
            [ s1( Select3 ), s2( Select3 ) ] = DDNormalize( s1( Select3 ), c3( Select3 ) );
            Select3A = Select3 & ( s2 ~= 0 );
            if any( Select3A, 'all' )
                [ s2( Select3A ), s3( Select3A ) ] = DDNormalize( s2( Select3A ), c4( Select3A ) );
            end
            Select3B = Select3 & ~( s2 ~= 0 );
            if any( Select3B, 'all' )
                [ s1( Select3B ), s2( Select3B ) ] = DDNormalize( s1( Select3B ), c4( Select3B ) );
            end
        end
    end

    Select4 = ~Select;
    if any( Select4, 'all' )
        [ s0( Select4 ), s1( Select4 ) ] = DDNormalize( s0( Select4 ), c2( Select4 ) );
        Select5 = Select4 & ( s1 ~= 0 );
        if any( Select5, 'all' )
            [ s1( Select5 ), s2( Select5 ) ] = DDNormalize( s1( Select5 ), c3( Select5 ) );
            Select5A = Select5 & ( s2 ~= 0 );
            if any( Select5A, 'all' )
                [ s2( Select5A ), s3( Select5A ) ] = DDNormalize( s2( Select5A ), c4( Select5A ) );
            end
            Select5B = Select5 & ~( s2 ~= 0 );
            if any( Select5B, 'all' )
                [ s1( Select5B ), s2( Select5B ) ] = DDNormalize( s1( Select5B ), c4( Select5B ) );
            end
        end
        Select6 = Select4 & ~( s1 ~= 0 );
        if any( Select6, 'all' )
            [ s0( Select6 ), s1( Select6 ) ] = DDNormalize( s0( Select6 ), c3( Select6 ) );
            Select6A = Select6 & ( s1 ~= 0 );
            if any( Select6A, 'all' )
                [ s1( Select6A ), s2( Select6A ) ] = DDNormalize( s1( Select6A ), c4( Select6A ) );
            end
            Select6B = Select6 & ~( s1 ~= 0 );
            if any( Select6B, 'all' )
                [ s0( Select6B ), s1( Select6B ) ] = DDNormalize( s0( Select6B ), c4( Select6B ) );
            end
        end
    end

    c0 = s0; c1 = s1; c2 = s2; c3 = s3;
end
