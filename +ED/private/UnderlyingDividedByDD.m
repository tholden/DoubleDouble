function [ r1, r2 ] = UnderlyingDividedByDD( a, b1, b2 )
    % Rescale to prevent overflow in intermediate products (cf. QD library)
    Rescale = abs( a ) > 2 ^ 969;
    if any( Rescale(:) )
        ScaleDown = 2 ^ -53;
        a( Rescale ) = a( Rescale ) * ScaleDown;
    end
    q1 = a ./ b1;
    [ p1, p2 ] = DDTimesUnderlying( b1, b2, q1 );
    [ r1, r2 ] = DDPlusUnderlying( -p1, -p2, a );
    q2 = r1 ./ b1;
    [ p1, p2 ] = DDTimesUnderlying( b1, b2, q2 );
    [ r1, ~  ] = DDPlusDD( r1, r2, -p1, -p2 );
    q3 = r1 ./ b1;
    [ q1, q2 ] = DDNormalize( q1, q2 );
    [ r1, r2 ] = DDPlusUnderlying( q1, q2, q3 );
    if any( Rescale(:) )
        ScaleUp = 2 ^ 53;
        r1( Rescale ) = r1( Rescale ) * ScaleUp;
        r2( Rescale ) = r2( Rescale ) * ScaleUp;
    end
    Select = ( b1 == 0 ) & ( b2 == 0 );
    if any( Select(:) )
        if ( isscalar( Select ) ) && ( numel( a ) > 1 )
            Select = repmat( Select, size( a ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( a ) )
            a = repmat( a, size( Select ) );
        end
        aSelect = a( Select );
        aSelect = sign( aSelect ) .* Inf;
        r1( Select ) = aSelect;
        r2( Select ) = aSelect;
    end
    Select = isinf( b1 );
    if any( Select(:) )
        if ( isscalar( Select ) ) && ( numel( a ) > 1 )
            Select = repmat( Select, size( a ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( a ) )
            a = repmat( a, size( Select ) );
        end
        aSelect = a( Select );
        aSelectSelect = ~isfinite( aSelect );
        aSelect = 0;
        aSelect( aSelectSelect ) = NaN;
        r1( Select ) = aSelect;
        r2( Select ) = aSelect;
    end
end
