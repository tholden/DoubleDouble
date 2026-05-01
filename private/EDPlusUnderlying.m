function [ s1, s2 ] = EDPlusUnderlying( a1, a2, b )
    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b );
    s2 = s2 + a2;
    [ s1, s2 ] = EDNormalize( s1, s2 );
end
