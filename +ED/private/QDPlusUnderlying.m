function [ c0, c1, c2, c3 ] = QDPlusUnderlying( a0, a1, a2, a3, b )
    [ c0, c1, c2, c3 ] = QDPlusQD( a0, a1, a2, a3, b, zeros( size( b ), 'like', b ), zeros( size( b ), 'like', b ), zeros( size( b ), 'like', b ) );
end
