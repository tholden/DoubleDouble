function [ p0, p1, p2, p3 ] = QDTimesDD( a0, a1, a2, a3, b0, b1 )
    [ p0, q0 ] = UnderlyingTimesUnderlying( a0, b0 );
    [ p1, q1 ] = UnderlyingTimesUnderlying( a0, b1 );
    [ p2, q2 ] = UnderlyingTimesUnderlying( a1, b0 );
    [ p3, q3 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p4, q4 ] = UnderlyingTimesUnderlying( a2, b0 );
    
    % Exact products for lower order terms
    [ p5, q5 ] = UnderlyingTimesUnderlying( a2, b1 );
    [ p6, q6 ] = UnderlyingTimesUnderlying( a3, b0 );
    [ p7, q7 ] = UnderlyingTimesUnderlying( a3, b1 );
    
    [ p1, p2, q0 ] = ThreeSum( p1, p2, q0 );
    [ p2, p3, p4 ] = ThreeSum( p2, p3, p4 );
    
    [ q1, q2, e_q12 ] = ThreeSum( q1, q2, 0 );
    
    [ s0, t0 ] = UnderlyingPlusUnderlying( p2, q1 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( p3, q2 );
    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    [ s2, e_s2 ] = UnderlyingPlusUnderlying( t0, t1 );
    [ s2, e_s2p4 ] = UnderlyingPlusUnderlying( s2, p4 );
    
    p2 = s0;
    
    [ p3, e1 ] = UnderlyingPlusUnderlying( p5, p6 );
    [ p3, e2 ] = UnderlyingPlusUnderlying( p3, q3 );
    [ p3, e3 ] = UnderlyingPlusUnderlying( p3, q4 );
    
    [ p3, q0, s1 ] = ThreeSum( p3, q0, s1 );
    [ p4_new, e_p4 ] = UnderlyingPlusUnderlying( q0, s2 );
    [ p3, p4_new ] = UnderlyingPlusUnderlying( p3, p4_new );
    
    [ p0, p1, p2, p3 ] = Renorm5( p0, p1, p2, p3, p4_new );
    
    % Cascade all tracked error terms accurately using a single QDPlusUnderlying
    % for performance, as they are all O(eps^4) or smaller.
    e_sum = e_s2 + e_s2p4 + e1 + e2 + e3 + e_p4 + e_q12 + q5 + q6 + p7 + q7;
    [ p0, p1, p2, p3 ] = QDPlusUnderlying( p0, p1, p2, p3, e_sum );
end
