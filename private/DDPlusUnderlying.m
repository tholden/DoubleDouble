function [ s1, s2 ] = DDPlusUnderlying( a1, a2, b )
    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b );
    s2 = s2 + a2;
    [ s1, s2 ] = DDNormalize( s1, s2 );
end
