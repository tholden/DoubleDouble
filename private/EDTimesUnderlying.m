function [ p1, p2 ] = EDTimesUnderlying( a1, a2, b )
    [ p1, p2 ] = UnderlyingTimesUnderlying( a1, b );
    p2 = p2 + a2 .* b;
    [ p1, p2 ] = EDNormalize( p1, p2 );
end
