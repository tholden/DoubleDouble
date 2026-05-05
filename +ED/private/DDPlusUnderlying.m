function [ s1, s2 ] = DDPlusUnderlying( a1, a2, b )
    [ s1, e1 ] = UnderlyingPlusUnderlying( a1, b );
    [ s2, e2 ] = UnderlyingPlusUnderlying( e1, a2 );
    [ s1, s2 ] = DDNormalize( s1, s2 );
    s2 = s2 + e2;
    [ s1, s2 ] = DDNormalize( s1, s2 );
end
