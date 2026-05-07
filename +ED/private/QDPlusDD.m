function [ s0, s1, s2, s3, s4 ] = QDPlusDD( a0, a1, a2, a3, b0, b1 )
    [ s0, t0 ] = UnderlyingPlusUnderlying( a0, b0 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    s2 = a2;
    [ s2, t0, t1 ] = ThreeSum( s2, t0, t1 );
    [ s3, t0 ] = UnderlyingPlusUnderlying( t0, a3 );
    [ t0, e1 ] = UnderlyingPlusUnderlying( t0, t1 );
    t0 = t0 + e1;
    [ s0, s1, s2, s3, ~ ] = Renorm5( s0, s1, s2, s3, t0 );
    s4 = 0;
end
