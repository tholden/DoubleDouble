function [ s0, s1, s2, s3, s4, s5, s6, s7 ] = QDTimesQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % Exact QD*QD multiplication producing 8 components.
    % Computes all 16 TwoProds (yielding 32 primitives) and exactly 
    % accumulates them via an unrolled Ripple-Insertion cascade.

    if ~isreal( a0 ) || ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( a3 ) || ...
            ~isreal( b0 ) || ~isreal( b1 ) || ~isreal( b2 ) || ~isreal( b3 )
        [ s0r, s1r, s2r, s3r, s4r, s5r, s6r, s7r ] = QDTimesQD( real( a0 ), real( a1 ), real( a2 ), real( a3 ),  ...
            real( b0 ), real( b1 ), real( b2 ), real( b3 ) );
        [ s0i, s1i, s2i, s3i, s4i, s5i, s6i, s7i ] = QDTimesQD( imag( a0 ), imag( a1 ), imag( a2 ), imag( a3 ),  ...
            imag( b0 ), imag( b1 ), imag( b2 ), imag( b3 ) );
        s0 = complex( s0r, s0i );
        s1 = complex( s1r, s1i );
        s2 = complex( s2r, s2i );
        s3 = complex( s3r, s3i );
        s4 = complex( s4r, s4i );
        s5 = complex( s5r, s5i );
        s6 = complex( s6r, s6i );
        s7 = complex( s7r, s7i );
        return
    end

    [ a0, a1, a2, a3, b0, b1, b2, b3 ] = ED.ExpandSingleton( a0, a1, a2, a3, b0, b1, b2, b3 );

    % O(eps^0)
    [ p00, q00 ] = UnderlyingTimesUnderlying( a0, b0 );

    % O(eps^1)
    [ p01, q01 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p10, q10 ] = UnderlyingTimesUnderlying( a1, b0 );

    % O(eps^2)
    [ p02, q02 ] = UnderlyingTimesUnderlying( a0, b2 );
    [ p11, q11 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p20, q20 ] = UnderlyingTimesUnderlying( a2, b0 );

    % O(eps^3)
    [ p03, q03 ] = UnderlyingTimesUnderlying( a0, b3 );
    [ p12, q12 ] = UnderlyingTimesUnderlying( a1, b2 );
    [ p21, q21 ] = UnderlyingTimesUnderlying( a2, b1 );
    [ p30, q30 ] = UnderlyingTimesUnderlying( a3, b0 );

    % O(eps^4)
    [ p13, q13 ] = UnderlyingTimesUnderlying( a1, b3 );
    [ p22, q22 ] = UnderlyingTimesUnderlying( a2, b2 );
    [ p31, q31 ] = UnderlyingTimesUnderlying( a3, b1 );

    % O(eps^5)
    [ p23, q23 ] = UnderlyingTimesUnderlying( a2, b3 );
    [ p32, q32 ] = UnderlyingTimesUnderlying( a3, b2 );

    % O(eps^6)
    [ p33, q33 ] = UnderlyingTimesUnderlying( a3, b3 );

    % Initialize 8-component state
    s0 = p00;
    s1 = q00;
    s2 = zeros( size( s0 ), 'like', s0 );
    s3 = s2; s4 = s2; s5 = s2; s6 = s2; s7 = s2;

    % List of 30 remaining terms, roughly ordered by magnitude
    terms = { p01, p10, q01, q10, ...
              p02, p11, p20, q02, q11, q20, ...
              p03, p12, p21, p30, q03, q12, q21, q30, ...
              p13, p22, p31, q13, q22, q31, ...
              p23, p32, q23, q32, ...
              p33, q33 };

    % Unrolled Ripple-Insertion Cascade using exact 6-op TwoSum
    for k = 1:30
        y = terms{k};
        
        s_new = s0 + y; bb = s_new - s0; y = (s0 - (s_new - bb)) + (y - bb); s0 = s_new;
        s_new = s1 + y; bb = s_new - s1; y = (s1 - (s_new - bb)) + (y - bb); s1 = s_new;
        s_new = s2 + y; bb = s_new - s2; y = (s2 - (s_new - bb)) + (y - bb); s2 = s_new;
        s_new = s3 + y; bb = s_new - s3; y = (s3 - (s_new - bb)) + (y - bb); s3 = s_new;
        s_new = s4 + y; bb = s_new - s4; y = (s4 - (s_new - bb)) + (y - bb); s4 = s_new;
        s_new = s5 + y; bb = s_new - s5; y = (s5 - (s_new - bb)) + (y - bb); s5 = s_new;
        s_new = s6 + y; bb = s_new - s6; y = (s6 - (s_new - bb)) + (y - bb); s6 = s_new;
        s_new = s7 + y; bb = s_new - s7; s7 = s_new; % Last remainder 'y' is discarded (falls beyond 8th component)
    end

    [ s0, s1, s2, s3 ] = QDNormalize( s0, s1, s2, s3 );
    [ s4, s5, s6, s7 ] = QDNormalize( s4, s5, s6, s7 );
    
    [ s3, s4 ] = DDNormalize( s3, s4 );
    [ s2, s3 ] = DDNormalize( s2, s3 );
    [ s1, s2 ] = DDNormalize( s1, s2 );
    [ s0, s1 ] = DDNormalize( s0, s1 );

end
