function [ p1, p2 ] = DDTimesUnderlying( a1, a2, b )
    % DD * underlying multiplication ( cf. operator*(dd_real, double) in QD library ).
    % Variable names match C++: p1,p2.
    % Improvement over C++: uses two_sum for p2 accumulation.
    [ p1, p2 ] = UnderlyingTimesUnderlying( a1, b );
    [ p2, e ] = UnderlyingPlusUnderlying( p2, a2 .* b );
    [ p1, p2 ] = DDNormalize( p1, p2 );
    p2 = p2 + e;
    [ p1, p2 ] = DDNormalize( p1, p2 );
end
