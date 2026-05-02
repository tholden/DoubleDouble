function [ s1, s2 ] = DDPlusDD( a1, a2, b1, b2 )
    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ t1, t2 ] = UnderlyingPlusUnderlying( a2, b2 );
    s2 = s2 + t1;
    [ s1, s2 ] = DDNormalize( s1, s2 );
    s2 = s2 + t2;
    [ s1, s2 ] = DDNormalize( s1, s2 );
end
