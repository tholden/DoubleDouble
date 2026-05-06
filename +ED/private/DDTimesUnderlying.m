function [ P1, P2, P3 ] = DDTimesUnderlying( A1, A2, A3, B )
    [ P1, P2 ] = UnderlyingTimesUnderlying( A1, B );
    [ P2, E1 ] = UnderlyingPlusUnderlying( P2, A2 .* B );
    [ P1, P2 ] = DDNormalize( P1, P2 );
    P2 = P2 + E1;
    [ P1, P2 ] = DDNormalize( P1, P2 );
    P3 = 0;
    if ~all( A3 == 0, 'all' )
        P3 = A3 .* B;
    end
end
