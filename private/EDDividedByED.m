function [ r1, r2 ] = EDDividedByED( a1, a2, b1, b2 )
    % Rescale to prevent overflow in intermediate products (cf. QD library)
    Rescale = abs( a1 ) > 2 ^ 969;
    if any( Rescale(:) )
        ScaleDown = 2 ^ -53;
        a1( Rescale ) = a1( Rescale ) * ScaleDown;
        a2( Rescale ) = a2( Rescale ) * ScaleDown;
    end
    q1 = a1 ./ b1;
    [ p1, p2 ] = EDTimesUnderlying( b1, b2, q1 );
    [ r1, r2 ] = EDPlusED( a1, a2, -p1, -p2 );
    q2 = r1 ./ b1;
    [ p1, p2 ] = EDTimesUnderlying( b1, b2, q2 );
    [ r1, ~  ] = EDPlusED( r1, r2, -p1, -p2 );
    q3 = r1 ./ b1;
    [ q1, q2 ] = EDNormalize( q1, q2 );
    [ r1, r2 ] = EDPlusED( q1, q2, q3, zeros( size( q3 ) ) );
    if any( Rescale(:) )
        ScaleUp = 2 ^ 53;
        r1( Rescale ) = r1( Rescale ) * ScaleUp;
        r2( Rescale ) = r2( Rescale ) * ScaleUp;
    end
    Select = ( b1 == 0 ) & ( b2 == 0 );
    if any( Select(:) )
        if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
            Select = repmat( Select, size( a1 ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
            a1 = repmat( a1, size( Select ) );
        end
        a1Select = a1( Select );
        a1Select = sign( a1Select ) .* Inf;
        r1( Select ) = a1Select;
        r2( Select ) = a1Select;
    end
    Select = isinf( b1 );
    if any( Select(:) )
        if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
            Select = repmat( Select, size( a1 ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
            a1 = repmat( a1, size( Select ) );
        end
        a1Select = a1( Select );
        a1SelectSelect = ~isfinite( a1Select );
        a1Select = 0;
        a1Select( a1SelectSelect ) = NaN;
        r1( Select ) = a1Select;
        r2( Select ) = a1Select;
    end
end
