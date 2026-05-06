function [ R1, R2, R3 ] = DDDividedByUnderlying( A1, A2, A3, B )
    % Rescale to prevent overflow in intermediate products ( cf. QD library )
    Rescale = abs( A1 ) > 2 ^ 969;
    if any( Rescale( : ) )
        ScaleDown = 2 ^ -53;
        A1( Rescale ) = A1( Rescale ) * ScaleDown;
        A2( Rescale ) = A2( Rescale ) * ScaleDown;
        A3( Rescale ) = A3( Rescale ) * ScaleDown;
    end
    R1 = A1 ./ B;
    [ P1, P2 ] = UnderlyingTimesUnderlying( R1, B );
    [ S, E ] = UnderlyingPlusUnderlying( A1, -P1 );
    E = E + A2;
    E = E - P2;
    T = S + E;
    R2 = T ./ B;
    R3 = A3 ./ B;
    [ R1, R2, R3 ] = DDNormalize( R1, R2, R3 );
    if any( Rescale( : ) )
        ScaleUp = 2 ^ 53;
        R1( Rescale ) = R1( Rescale ) * ScaleUp;
        R2( Rescale ) = R2( Rescale ) * ScaleUp;
        R3( Rescale ) = R3( Rescale ) * ScaleUp;
    end
    Select = B == 0;
    if any( Select( : ) )
        if ( isscalar( Select ) ) && ( numel( A1 ) > 1 )
            Select = repmat( Select, size( A1 ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( A1 ) )
            A1 = repmat( A1, size( Select ) );
        end
        A1Select = A1( Select );
        A1Select = sign( A1Select ) .* Inf;
        R1( Select ) = A1Select;
        R2( Select ) = A1Select;
        R3( Select ) = 0;
    end
    Select = isinf( B );
    if any( Select( : ) )
        if ( isscalar( Select ) ) && ( numel( A1 ) > 1 )
            Select = repmat( Select, size( A1 ) );
        elseif ( numel( Select ) > 1 ) && ( isscalar( A1 ) )
            A1 = repmat( A1, size( Select ) );
        end
        A1Select = A1( Select );
        A1SelectSelect = ~isfinite( A1Select );
        A1Select = 0;
        A1Select( A1SelectSelect ) = NaN;
        R1( Select ) = A1Select;
        R2( Select ) = A1Select;
        R3( Select ) = 0;
    end
end
