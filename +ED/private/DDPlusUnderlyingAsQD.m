function [ c0, c1, c2, c3 ] = DDPlusUnderlyingAsQD( a0, a1, b )
    [ c0, e ] = UnderlyingPlusUnderlying( a0, b );
    [ c1, c2 ] = UnderlyingPlusUnderlying( a1, e );
    c3 = 0;
    [ c0, c1, c2, c3 ] = QDNormalize( c0, c1, c2, c3 );
end
