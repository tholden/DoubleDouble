function [ p1, p2 ] = DDTimesDD( a1, a2, b1, b2 )
    [ p1, p2 ] = UnderlyingTimesUnderlying( a1, b1 );
    t = a1 .* b2 + a2 .* b1;
    p2 = p2 + t;
    [ p1, p2 ] = DDNormalize( p1, p2 );
end
