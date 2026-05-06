function [ r1, r2 ] = DDDividedByUnderlying( a1, a2, b )
    % Rescale to prevent overflow in intermediate products ( cf. QD library )
    Rescale = abs( a1 ) > 2 ^ 969;
    if any( Rescale( : ) )
        ScaleDown = 2 ^ -53;
        a1( Rescale ) = a1( Rescale ) * ScaleDown;
        a2( Rescale ) = a2( Rescale ) * ScaleDown;
    end
    r1 = a1 ./ b;
    [ p1, p2 ] = UnderlyingTimesUnderlying( r1, b );
    [ s, e ] = UnderlyingPlusUnderlying( a1, -p1 );
    e = e + a2;
    e = e - p2;
    t = s + e;
    r2 = t ./ b;
    [ r1, r2 ] = DDNormalize( r1, r2 );
    if any( Rescale( : ) )
        ScaleUp = 2 ^ 53;
        r1( Rescale ) = r1( Rescale ) * ScaleUp;
        r2( Rescale ) = r2( Rescale ) * ScaleUp;
    end
    Select = b == 0;
    if any( Select( : ) )
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
    Select = isinf( b );
    if any( Select( : ) )
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
