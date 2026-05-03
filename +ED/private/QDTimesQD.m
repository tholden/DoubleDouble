function [ s0, s1, s2, s3 ] = QDTimesQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % Accurate multiplication ( cf. qd_real::accurate_mul in QD library ).
    % Uses 10 UnderlyingTimesUnderlyings and Nine-Two-Sum for O( eps^3 ) accumulation.
    % O( eps^1 ) and O( eps^2 ) terms
    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b0 );
    [ p1, q1 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a1, b0 );
    [ p3, q3 ] = UnderlyingTimesUnderlying( a0, b2 );
    [ p4, q4 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p5, q5 ] = UnderlyingTimesUnderlying( a2, b0 );
    % Start accumulation
    [ p1, p2, q0 ] = ThreeSum( p1, p2, q0 );
    % Six-Three-Sum of p2, q1, q2, p3, p4, p5
    [ p2, q1, q2 ] = ThreeSum( p2, q1, q2 );
    [ p3, p4, p5 ] = ThreeSum( p3, p4, p5 );
    [ ss0, t0 ] = UnderlyingPlusUnderlying( p2, p3 );
    [ ss1, t1 ] = UnderlyingPlusUnderlying( q1, p4 );
    ss2 = q2 + p5;
    [ ss1, t0_2 ] = UnderlyingPlusUnderlying( ss1, t0 );
    ss2 = ss2 + ( t0_2 + t1 );
    % O( eps^3 ) terms
    [ p6, q6 ] = UnderlyingTimesUnderlying( a0, b3 );
    [ p7, q7 ] = UnderlyingTimesUnderlying( a1, b2 );
    [ p8, q8 ] = UnderlyingTimesUnderlying( a2, b1 );
    [ p9, q9 ] = UnderlyingTimesUnderlying( a3, b0 );
    % Nine-Two-Sum of q0, ss1, q3, q4, q5, p6, p7, p8, p9
    [ q0, q3 ] = UnderlyingPlusUnderlying( q0, q3 );
    [ q4, q5 ] = UnderlyingPlusUnderlying( q4, q5 );
    [ p6, p7 ] = UnderlyingPlusUnderlying( p6, p7 );
    [ p8, p9 ] = UnderlyingPlusUnderlying( p8, p9 );
    [ t0, t1 ] = UnderlyingPlusUnderlying( q0, q4 );
    t1 = t1 + ( q3 + q5 );
    [ r0, r1 ] = UnderlyingPlusUnderlying( p6, p8 );
    r1 = r1 + ( p7 + p9 );
    [ q3_2, q4_2 ] = UnderlyingPlusUnderlying( t0, r0 );
    q4_2 = q4_2 + ( t1 + r1 );
    [ t0, t1 ] = UnderlyingPlusUnderlying( q3_2, ss1 );
    t1 = t1 + q4_2;
    % O( eps^4 ) terms -- Nine-One-Sum
    t1 = t1 + a1.*b3 + a2.*b2 + a3.*b1 + q6 + q7 + q8 + q9 + ss2;
    [ s0, s1, s2, s3 ] = Renorm5( p0, p1, ss0, t0, t1 );
end
