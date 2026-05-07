function [ s0, s1, s2, s3 ] = QDTimesUnderlying( a0, a1, a2, a3, b )
    % QD * underlying ( cf. operator*(qd_real, double) in QD library ).

    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b );
    [ p1, q1 ] = UnderlyingTimesUnderlying( a1, b );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a2, b );
    p3 = a3 .* b;

    s0 = p0;

    [ s1, s2 ] = UnderlyingPlusUnderlying( q0, p1 );

    [ s2, q1, p2 ] = ThreeSum( s2, q1, p2 );

    [ q1, q2 ] = ThreeSum2( q1, q2, p3 );
    s3 = q1;

    s4 = q2 + p2;

    [ s0, s1, s2, s3 ] = Renorm5( s0, s1, s2, s3, s4 );
end
