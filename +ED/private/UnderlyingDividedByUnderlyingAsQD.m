function [ s0, s1, s2, s3 ] = UnderlyingDividedByUnderlyingAsQD( a, b )
    Rescale = abs( a ) > 2 ^ 969;
    if any( Rescale, 'all' )
        ScaleDown = 2 ^ -53;
        a( Rescale ) = a( Rescale ) * ScaleDown;
    end
    q0 = a ./ b;
    [ r0, r1 ] = UnderlyingTimesUnderlying( b, q0 );
    [ r0, r1, r2, r3 ] = DDPlusUnderlyingAsQD( -r0, -r1, a );
    q1 = r0 ./ b;
    [ t0, t1 ] = UnderlyingTimesUnderlying( b, q1 );
    [ r0, r1, r2, r3 ] = QDPlusDD( r0, r1, r2, r3, -t0, -t1 );
    q2 = r0 ./ b;
    [ t0, t1 ] = UnderlyingTimesUnderlying( b, q2 );
    [ r0, r1, r2, r3 ] = QDPlusDD( r0, r1, r2, r3, -t0, -t1 );
    q3 = r0 ./ b;
    [ t0, t1 ] = UnderlyingTimesUnderlying( b, q3 );
    [ r0, ~, ~, ~ ] = QDPlusDD( r0, r1, r2, r3, -t0, -t1 );
    q4 = r0 ./ b;
    [ s0, s1, s2, s3 ] = Renorm5( q0, q1, q2, q3, q4 );
    if any( Rescale, 'all' )
        ScaleUp = 2 ^ 53;
        s0( Rescale ) = s0( Rescale ) * ScaleUp;
        s1( Rescale ) = s1( Rescale ) * ScaleUp;
        s2( Rescale ) = s2( Rescale ) * ScaleUp;
        s3( Rescale ) = s3( Rescale ) * ScaleUp;
    end
end
