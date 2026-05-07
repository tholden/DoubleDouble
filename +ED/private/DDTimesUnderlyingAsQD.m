function [ p0, p1, p2, p3, p4 ] = DDTimesUnderlyingAsQD( a0, a1, b0 )
    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b0 );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a1, b0 );
    [ p1, p2 ] = UnderlyingPlusUnderlying( p2, q0 );
    [ p2, p3 ] = UnderlyingPlusUnderlying( p2, q2 );
    [ p0, p1, p2, p3, p4 ] = QDNormalize( p0, p1, p2, p3 );
end
