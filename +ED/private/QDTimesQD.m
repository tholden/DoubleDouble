function [ s0, s1, s2, s3, s4 ] = QDTimesQD( a0, a1, a2, a3, b0, b1, b2, b3, a4, b4 )
    % Accurate multiplication ( cf. qd_real::accurate_mul in QD library ).
    % a4, b4 are optional 5th components (v3) from BaseQuadDouble.
    if nargin < 9
        a4 = 0;
    end
    if nargin < 10
        b4 = 0;
    end
    
    % Calculate all 13 possible exact products up to O(eps^4)
    [ p00, q00 ] = UnderlyingTimesUnderlying( a0, b0 );
    [ p01, q01 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p10, q10 ] = UnderlyingTimesUnderlying( a1, b0 );
    
    [ p02, q02 ] = UnderlyingTimesUnderlying( a0, b2 );
    [ p11, q11 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p20, q20 ] = UnderlyingTimesUnderlying( a2, b0 );
    
    [ p03, q03 ] = UnderlyingTimesUnderlying( a0, b3 );
    [ p12, q12 ] = UnderlyingTimesUnderlying( a1, b2 );
    [ p21, q21 ] = UnderlyingTimesUnderlying( a2, b1 );
    [ p30, q30 ] = UnderlyingTimesUnderlying( a3, b0 );
    
    [ p13, q13 ] = UnderlyingTimesUnderlying( a1, b3 );
    [ p22, q22 ] = UnderlyingTimesUnderlying( a2, b2 );
    [ p31, q31 ] = UnderlyingTimesUnderlying( a3, b1 );

    % Use the fast cascade for the dominant terms to establish the base QD
    [ p01, p10, q00 ] = ThreeSum( p01, p10, q00 );
    
    [ p10, q01, q10 ] = ThreeSum( p10, q01, q10 );
    [ p02, p11, p20 ] = ThreeSum( p02, p11, p20 );
    
    [ ss0, t0 ] = UnderlyingPlusUnderlying( p10, p02 );
    [ ss1, t1 ] = UnderlyingPlusUnderlying( q01, p11 );
    [ ss2, e_ss2 ] = UnderlyingPlusUnderlying( q10, p20 );
    
    [ ss1, t0_2 ] = UnderlyingPlusUnderlying( ss1, t0 );
    [ ss2_t, e_t ] = UnderlyingPlusUnderlying( t0_2, t1 );
    [ ss2, e_ss2b ] = UnderlyingPlusUnderlying( ss2, ss2_t );
    
    [ q00, p03 ] = UnderlyingPlusUnderlying( q00, p03 );
    [ p12, p21 ] = UnderlyingPlusUnderlying( p12, p21 );
    [ p30, p13 ] = UnderlyingPlusUnderlying( p30, p13 );
    [ p22, p31 ] = UnderlyingPlusUnderlying( p22, p31 );
    
    [ t0, t1 ] = UnderlyingPlusUnderlying( q00, p12 );
    [ t1_sum, e_t1 ] = UnderlyingPlusUnderlying( t1, p03 );
    [ t1_sum, e_t1b ] = UnderlyingPlusUnderlying( t1_sum, p21 );
    
    [ r0, r1 ] = UnderlyingPlusUnderlying( p30, p22 );
    [ r1_sum, e_r1 ] = UnderlyingPlusUnderlying( r1, p13 );
    [ r1_sum, e_r1b ] = UnderlyingPlusUnderlying( r1_sum, p31 );
    
    [ q3_2, q4_2 ] = UnderlyingPlusUnderlying( t0, r0 );
    [ q4_sum, e_q4 ] = UnderlyingPlusUnderlying( q4_2, t1_sum );
    [ q4_sum, e_q4b ] = UnderlyingPlusUnderlying( q4_sum, r1_sum );
    
    [ t0, t1 ] = UnderlyingPlusUnderlying( q3_2, ss1 );
    [ t1, e_t1c ] = UnderlyingPlusUnderlying( t1, q4_sum );
    
    % Base QuadDouble
    [ p0, p1, p2, p3, ~ ] = Renorm5( p00, p01, ss0, t0, t1 );
    
    % ss2 is O(eps^2), so it must be added separately to prevent swallowing 
    % the O(eps^4) terms during sloppy addition.
    [ p0, p1, p2, p3, ~ ] = QDPlusUnderlying( p0, p1, p2, p3, ss2 );
    
    % Cascade all remaining O(eps^3) and smaller tracked error terms.
    % Include v3 cross-terms: a4*b0 + a0*b4 contribute at O(eps^4).
    e_sum = e_ss2 + e_ss2b + e_t + e_t1 + e_t1b + e_r1 + e_r1b + e_q4 + e_q4b + e_t1c + ...
            q02 + q11 + q20 + q03 + q12 + q21 + q30 + q13 + q22 + q31 + ...
            a4 .* b0 + a0 .* b4;
            
    [ s0, s1, s2, s3, ~ ] = QDPlusUnderlying( p0, p1, p2, p3, e_sum );
    s4 = 0;
end
