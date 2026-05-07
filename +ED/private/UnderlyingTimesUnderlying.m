function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
    p1 = a .* b;
    [ a1, a2 ] = Split( a );
    [ b1, b2 ] = Split( b );
    [ t1, e1 ] = UnderlyingPlusUnderlying( a1 .* b1, -p1 );
    [ t2, e2 ] = UnderlyingPlusUnderlying( t1, a1 .* b2 );
    [ t3, e3 ] = UnderlyingPlusUnderlying( t2, a2 .* b1 );
    p2 = t3 + (a2 .* b2 + e1 + e2 + e3);
end
