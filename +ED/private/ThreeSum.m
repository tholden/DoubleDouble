function [ a, b, c ] = ThreeSum( a, b, c )
    [ t1, t2 ] = UnderlyingPlusUnderlying( a, b );
    [ a, t3 ]  = UnderlyingPlusUnderlying( c, t1 );
    [ b, c ]   = UnderlyingPlusUnderlying( t2, t3 );
end
