function [ P1, P2, P3 ] = DDTimesDD( A1, A2, A3, B1, B2, B3 )
    [ P1, P2 ] = UnderlyingTimesUnderlying( A1, B1 );
    [ C1, C2 ] = UnderlyingPlusUnderlying( P2, A1 .* B2 );
    [ P2, E1 ] = UnderlyingPlusUnderlying( C1, A2 .* B1 );
    [ P1, P2 ] = DDNormalize( P1, P2 );
    Correction = C2 + E1 + A2 .* B2;
    P2Old = P2;
    P2 = P2 + Correction;
    [ P1, P2 ] = DDNormalize( P1, P2 );
    P3 = Correction - ( P2 - P2Old );
    if ~all( A3 == 0, 'all' )
        P3 = P3 + A1 .* B3;
    end
    if ~all( B3 == 0, 'all' )
        P3 = P3 + A3 .* B1;
    end
end
