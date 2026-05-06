function [ S1, S2, S3 ] = DDPlusUnderlying( A1, A2, A3, B )
    [ S1, E1 ] = UnderlyingPlusUnderlying( A1, B );
    [ S2, E2 ] = UnderlyingPlusUnderlying( E1, A2 );
    [ S1, S2 ] = DDNormalize( S1, S2 );
    S2 = S2 + E2;
    [ S1, S2 ] = DDNormalize( S1, S2 );
    S3 = A3;
end
