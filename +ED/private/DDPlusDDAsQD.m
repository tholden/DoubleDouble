function [ s0, s1, s2, s3 ] = DDPlusDDAsQD( a0, a1, b0, b1 )
    [ s0, t0 ] = UnderlyingPlusUnderlying( a0, b0 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    [ s2, s3 ] = UnderlyingPlusUnderlying( t0, t1 );
    [ s0, s1, s2, s3 ] = QDNormalize( s0, s1, s2, s3 );
end
