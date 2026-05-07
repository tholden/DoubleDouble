function [ c0, c1, c2, c3 ] = QDPlusUnderlying( a0, a1, a2, a3, b )
    [ c0, e ] = UnderlyingPlusUnderlying( a0, b );
    [ c1, e ] = UnderlyingPlusUnderlying( a1, e );
    [ c2, e ] = UnderlyingPlusUnderlying( a2, e );
    [ c3, e ] = UnderlyingPlusUnderlying( a3, e );
    [ c0, c1, c2, c3 ] = Renorm5( c0, c1, c2, c3, e );
end
