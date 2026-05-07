function [ p1, p2 ] = DDTimesDD( a1, a2, b1, b2 )
    % DD*DD multiplication ( cf. operator*(dd_real, dd_real) in QD library ).
    % Variable names match C++: p1,p2.
    % Improvement over C++: uses two_sum for intermediate accumulation.
    [ p1, p2 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ p2, e1 ] = UnderlyingPlusUnderlying( p2, a1 .* b2 );
    [ p2, e2 ] = UnderlyingPlusUnderlying( p2, a2 .* b1 );
    [ p1, p2 ] = DDNormalize( p1, p2 );
    p2 = p2 + (e1 + e2 + a2 .* b2);
    [ p1, p2 ] = DDNormalize( p1, p2 );
end
