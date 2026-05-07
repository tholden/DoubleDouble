function [ s0, s1, s2, s3 ] = QDPlusDD( a0, a1, a2, a3, b0, b1 )
    [ s0, s1, s2, s3 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, zeros( size( b1 ), 'like', b1 ), zeros( size( b1 ), 'like', b1 ) );
end
