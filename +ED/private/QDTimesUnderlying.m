function [ s0, s1, s2, s3, s4 ] = QDTimesUnderlying( a0, a1, a2, a3, b )
    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b );
    [ p1, q1 ] = UnderlyingTimesUnderlying( a1, b );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a2, b );
    p3 = a3 .* b;
    ss0 = p0;
    [ ss1, s2_1 ] = UnderlyingPlusUnderlying( q0, p1 );
    [ s2_1, q1, p2 ] = ThreeSum( s2_1, q1, p2 );
    [ q1, q2, p3_err ] = ThreeSum( q1, q2, p3 );
    ss3 = q1;
    [ s4, e_s4 ] = UnderlyingPlusUnderlying( q2, p2 );
    s4 = s4 + ( e_s4 + p3_err );
    [ s0, s1, s2, s3, s4 ] = Renorm5( ss0, ss1, s2_1, ss3, s4 );
end
