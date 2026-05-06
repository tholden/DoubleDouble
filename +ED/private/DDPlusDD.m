function [ S1, S2, S3 ] = DDPlusDD( A1, A2, A3, B1, B2, B3 )
    [ S1, E1 ] = UnderlyingPlusUnderlying( A1, B1 );
    [ T1, E2 ] = UnderlyingPlusUnderlying( A2, B2 );
    [ S2, E3 ] = UnderlyingPlusUnderlying( E1, T1 );
    [ S1, S2 ] = DDNormalize( S1, S2 );
    S2 = S2 + ( E2 + E3 );
    [ S1, S2 ] = DDNormalize( S1, S2 );
    S3 = A3 + B3;
end
