function [ s1, s2 ] = DDNormalize( a1, a2 )
    if isobject( a1 )
        [ s1, s2 ] = a1.UnderlyingPlusUnderlying( a2 );
    else
        a1 = Normalize( a1 );
        a2 = Normalize( a2 );
        s1 = a1 + a2;
        t = s1 - a1;
        s2 = a2 - t;
    end
end
