function [ a, b ] = ThreeSum2( a, b, c )
    [ t1, t2 ] = UnderlyingPlusUnderlying( a, b );
    [ a, t3 ]  = UnderlyingPlusUnderlying( c, t1 );
    b = t2 + t3;
end
