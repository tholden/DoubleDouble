function [ p1, p2, p3, p4 ] = UnderlyingTimesUnderlyingAsQD( a, b )
    [ p1, p2 ] = UnderlyingTimesUnderlying( a, b );
    p3 = zeros( size( p1 ), 'like', p1 );
    p4 = zeros( size( p1 ), 'like', p1 );
end
