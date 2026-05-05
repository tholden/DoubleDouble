function [ p1, p2 ] = DDTimesDD( a1, a2, b1, b2 )
    [ p1, p2 ] = UnderlyingTimesUnderlying( a1, b1 );
    [ c1, c2 ] = UnderlyingPlusUnderlying( p2, a1 .* b2 );
    [ p2, e1 ] = UnderlyingPlusUnderlying( c1, a2 .* b1 );
    [ p1, p2 ] = DDNormalize( p1, p2 );
    p2 = p2 + (c2 + e1 + a2 .* b2);
    [ p1, p2 ] = DDNormalize( p1, p2 );
end
