function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
    p1 = a .* b;
    [ a1, a2 ] = Split( a );
    [ b1, b2 ] = Split( b );
    t1 = a1 .* b1 - p1;
    t2 = t1 + a1 .* b2 + a2 .* b1;
    p2 = t2 + a2 .* b2;
end
