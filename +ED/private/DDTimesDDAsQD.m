function [ p0, p1, p2, p3 ] = DDTimesDDAsQD( a0, a1, b0, b1 )
    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b0 );
    [ p1, q1 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a1, b0 );
    [ p3, q3 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p1, p2, q0 ] = ThreeSum( p1, p2, q0 );
    [ p2, p3 ] = UnderlyingPlusUnderlying( p2, p3 );
    [ q1, q2 ] = UnderlyingPlusUnderlying( q1, q2 );
    [ s0, t0 ] = UnderlyingPlusUnderlying( p2, q1 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( p3, q2 );
    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    s2 = t0 + t1;
    p2 = s0;
    [ p3, q0 ] = ThreeSum2( q3, q0, s1 );
    p4 = q0 + s2;
    [ p0, p1, p2, p3 ] = Renorm5( p0, p1, p2, p3, p4 );
end
