function [ s0, s1, s2, s3, s4, s5, s6, s7 ] = ODTimesOD( a0, a1, a2, a3, a4, a5, a6, a7, b0, b1, b2, b3, b4, b5, b6, b7 )
    % OctoDouble multiplication on 8+8 raw doubles.
    % Extends QDTimesQD (accurate_mul) from 4-slot to 8-slot precision.
    %
    % Slots 0-3 are computed identically to QDTimesQD.
    % Slots 4-7 extend the pattern with additional TwoProd terms.

    % ====== Identical to QDTimesQD: compute slots 0-3 ======

    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b0 );

    [ p1, q1 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a1, b0 );

    [ p3, q3 ] = UnderlyingTimesUnderlying( a0, b2 );
    [ p4, q4 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p5, q5 ] = UnderlyingTimesUnderlying( a2, b0 );

    % Start accumulation
    [ p1, p2, q0 ] = ThreeSum( p1, p2, q0 );

    % Six-Three-Sum of p2, q1, q2, p3, p4, p5.
    [ p2, q1, q2 ] = ThreeSum( p2, q1, q2 );
    [ p3, p4, p5 ] = ThreeSum( p3, p4, p5 );
    [ acc2, t0 ] = UnderlyingPlusUnderlying( p2, p3 );
    [ acc2b, t1 ] = UnderlyingPlusUnderlying( q1, p4 );
    [ acc2c, acc2d ] = UnderlyingPlusUnderlying( q2, p5 );
    [ acc2b, t0 ] = UnderlyingPlusUnderlying( acc2b, t0 );
    [ acc2c, acc2e ] = UnderlyingPlusUnderlying( acc2c, t0 );
    acc2d = acc2d + acc2e + t1;

    % O(eps^3) order terms
    [ p6, q6 ] = UnderlyingTimesUnderlying( a0, b3 );
    [ p7, q7 ] = UnderlyingTimesUnderlying( a1, b2 );
    [ p8, q8 ] = UnderlyingTimesUnderlying( a2, b1 );
    [ p9, q9 ] = UnderlyingTimesUnderlying( a3, b0 );

    % Nine-Two-Sum of q0, acc2b, q3, q4, q5, p6, p7, p8, p9.
    [ q0, q3 ] = UnderlyingPlusUnderlying( q0, q3 );
    [ q4, q5 ] = UnderlyingPlusUnderlying( q4, q5 );
    [ p6, p7 ] = UnderlyingPlusUnderlying( p6, p7 );
    [ p8, p9 ] = UnderlyingPlusUnderlying( p8, p9 );
    [ t0, t1 ] = UnderlyingPlusUnderlying( q0, q4 );
    [ t1, t1e ] = UnderlyingPlusUnderlying( t1, q3 );
    t1e = t1e + q5;
    [ r0, r1 ] = UnderlyingPlusUnderlying( p6, p8 );
    [ r1, r1e ] = UnderlyingPlusUnderlying( r1, p7 );
    r1e = r1e + p9;
    [ acc3_a, acc3_b ] = UnderlyingPlusUnderlying( t0, r0 );
    [ acc3_b, acc3_be ] = UnderlyingPlusUnderlying( acc3_b, t1 );
    acc3_be = acc3_be + r1 + t1e + r1e;
    [ acc3, acc3_err ] = UnderlyingPlusUnderlying( acc3_a, acc2b );
    [ acc3_err, acc3_err2 ] = UnderlyingPlusUnderlying( acc3_err, acc3_b );
    acc3_err2 = acc3_err2 + acc3_be;

    % ====== Slot 4: O(eps^4) terms -- 5 TwoProds ======
    [ pa0, qa0 ] = UnderlyingTimesUnderlying( a0, b4 );
    [ pa1, qa1 ] = UnderlyingTimesUnderlying( a1, b3 );
    [ pa2, qa2 ] = UnderlyingTimesUnderlying( a2, b2 );
    [ pa3, qa3 ] = UnderlyingTimesUnderlying( a3, b1 );
    [ pa4, qa4 ] = UnderlyingTimesUnderlying( a4, b0 );

    % Accumulate slot 4 with exact pairwise TwoSum pattern
    [ pa0, pa1 ] = UnderlyingPlusUnderlying( pa0, pa1 );
    [ pa2, pa3 ] = UnderlyingPlusUnderlying( pa2, pa3 );
    [ q6, q7 ] = UnderlyingPlusUnderlying( q6, q7 );
    [ q8, q9 ] = UnderlyingPlusUnderlying( q8, q9 );

    [ u0, u1 ] = UnderlyingPlusUnderlying( pa0, pa2 );
    [ u1, u1e ] = UnderlyingPlusUnderlying( u1, pa1 );
    u1e = u1e + pa3 + pa4;
    [ v0, v1 ] = UnderlyingPlusUnderlying( q6, q8 );
    [ v1, v1e ] = UnderlyingPlusUnderlying( v1, q7 );
    v1e = v1e + q9;
    [ w0, w1 ] = UnderlyingPlusUnderlying( u0, v0 );
    [ w1, w1e ] = UnderlyingPlusUnderlying( w1, u1 );
    w1e = w1e + v1 + u1e + v1e;

    [ acc4, acc4_tmp ] = UnderlyingPlusUnderlying( w0, acc3_err );
    [ acc4_tmp2, acc4_tmp3 ] = UnderlyingPlusUnderlying( acc4_tmp, acc2c );
    [ acc4_err, acc4_err2 ] = UnderlyingPlusUnderlying( acc4_tmp2, w1 );
    acc4_err = acc4_err + acc4_err2 + acc4_tmp3 + w1e + acc3_err2 + acc2d;

    % ====== Slot 5: O(eps^5) terms -- 6 TwoProds ======
    [ pb0, qb0 ] = UnderlyingTimesUnderlying( a0, b5 );
    [ pb1, qb1 ] = UnderlyingTimesUnderlying( a1, b4 );
    [ pb2, qb2 ] = UnderlyingTimesUnderlying( a2, b3 );
    [ pb3, qb3 ] = UnderlyingTimesUnderlying( a3, b2 );
    [ pb4, qb4 ] = UnderlyingTimesUnderlying( a4, b1 );
    [ pb5, qb5 ] = UnderlyingTimesUnderlying( a5, b0 );

    % Accumulate slot 5 with pairwise TwoSum pattern
    [ pb0, pb1 ] = UnderlyingPlusUnderlying( pb0, pb1 );
    [ pb2, pb3 ] = UnderlyingPlusUnderlying( pb2, pb3 );
    [ pb4, pb5 ] = UnderlyingPlusUnderlying( pb4, pb5 );
    [ qa0, qa1 ] = UnderlyingPlusUnderlying( qa0, qa1 );
    [ qa2, qa3 ] = UnderlyingPlusUnderlying( qa2, qa3 );

    [ u0, u1 ] = UnderlyingPlusUnderlying( pb0, pb2 );
    tmp = pb1 + pb3;
    u1 = u1 + tmp;
    [ v0, v1 ] = UnderlyingPlusUnderlying( qa0, qa2 );
    tmp = qa1 + qa3 + qa4;
    v1 = v1 + tmp;
    [ w0, w1 ] = UnderlyingPlusUnderlying( u0, v0 );
    tmp = u1 + v1;
    w1 = w1 + tmp;
    [ x0, x1 ] = UnderlyingPlusUnderlying( w0, pb4 );
    tmp = x1 + w1 + pb5;

    [ acc5, acc5_err ] = UnderlyingPlusUnderlying( x0, acc4_err );
    acc5_err = acc5_err + tmp;

    % ====== Slot 6: O(eps^6) -- sloppy ======
    tmp = a0 .* b6 + a1 .* b5 + a2 .* b4 + a3 .* b3 + a4 .* b2 + a5 .* b1 + a6 .* b0 ...
        + qb0 + qb1 + qb2 + qb3 + qb4 + qb5;
    acc6 = tmp + acc5_err;

    % ====== Slot 7: O(eps^7) -- sloppy ======
    acc7 = a0 .* b7 + a1 .* b6 + a2 .* b5 + a3 .* b4 ...
         + a4 .* b3 + a5 .* b2 + a6 .* b1 + a7 .* b0;

    % ====== Renormalize all 8 slots ======
    % Three-pass renormalization with full two_sum.
    % Pass 1: upward with full two_sum to fix magnitude ordering
    [ acc6, acc7 ] = UnderlyingPlusUnderlying( acc6, acc7 );
    [ acc5, acc6 ] = UnderlyingPlusUnderlying( acc5, acc6 );
    [ acc4, acc5 ] = UnderlyingPlusUnderlying( acc4, acc5 );
    [ acc3, acc4 ] = UnderlyingPlusUnderlying( acc3, acc4 );
    [ acc2, acc3 ] = UnderlyingPlusUnderlying( acc2, acc3 );
    [ p1, acc2 ]   = UnderlyingPlusUnderlying( p1, acc2 );
    [ p0, p1 ]     = UnderlyingPlusUnderlying( p0, p1 );
    % Pass 2: downward
    [ p0, p1 ]     = DDNormalize( p0, p1 );
    [ p1, acc2 ]   = DDNormalize( p1, acc2 );
    [ acc2, acc3 ] = DDNormalize( acc2, acc3 );
    [ acc3, acc4 ] = DDNormalize( acc3, acc4 );
    [ acc4, acc5 ] = DDNormalize( acc4, acc5 );
    [ acc5, acc6 ] = DDNormalize( acc5, acc6 );
    [ acc6, acc7 ] = DDNormalize( acc6, acc7 );
    % Pass 3: upward again
    [ acc6, acc7 ] = UnderlyingPlusUnderlying( acc6, acc7 );
    [ acc5, acc6 ] = UnderlyingPlusUnderlying( acc5, acc6 );
    [ acc4, acc5 ] = UnderlyingPlusUnderlying( acc4, acc5 );
    [ acc3, acc4 ] = UnderlyingPlusUnderlying( acc3, acc4 );
    [ acc2, acc3 ] = UnderlyingPlusUnderlying( acc2, acc3 );
    [ p1, acc2 ]   = UnderlyingPlusUnderlying( p1, acc2 );
    [ s0, s1 ]     = UnderlyingPlusUnderlying( p0, p1 );
    s2 = acc2; s3 = acc3; s4 = acc4; s5 = acc5; s6 = acc6; s7 = acc7;
end
