function [ R1, R2, R3 ] = UnderlyingDividedByDD( A1, B1, B2, B3 )
    % Rescale to prevent overflow in intermediate products ( cf. QD library )
    Rescale = abs( A1 ) > 2 ^ 969;
    if any( Rescale( : ) )
        ScaleDown = 2 ^ -53;
        A1( Rescale ) = A1( Rescale ) * ScaleDown;
    end
    Q1 = A1 ./ B1;
    [ P1, P2, ~ ] = DDTimesUnderlying( B1, B2, B3, Q1 );
    [ R1, R2, ~ ] = DDPlusUnderlying( -P1, -P2, 0, A1 );
    Q2 = R1 ./ B1;
    [ P1, P2, ~ ] = DDTimesUnderlying( B1, B2, B3, Q2 );
    [ R1, ~, ~ ] = DDPlusDD( R1, R2, 0, -P1, -P2, 0 );
    Q3 = R1 ./ B1;
    [ Q1, Q2 ] = DDNormalize( Q1, Q2 );
    [ R1, R2, R3 ] = DDPlusUnderlying( Q1, Q2, 0, Q3 );
    if any( Rescale( : ) )
        ScaleUp = 2 ^ 53;
        R1( Rescale ) = R1( Rescale ) * ScaleUp;
        R2( Rescale ) = R2( Rescale ) * ScaleUp;
        R3( Rescale ) = R3( Rescale ) * ScaleUp;
    end
    Select = ( B1 == 0 ) & ( B2 == 0 );
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
    Select = isinf( B1 );
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
