function [ s0, s1, s2, s3 ] = QDDividedByDD( a0, a1, a2, a3, b0, b1 )
    Rescale = abs( a0 ) > 2 ^ 969;
    if any( Rescale, 'all' )
        ScaleDown = 2 ^ -53;
        a0( Rescale ) = a0( Rescale ) * ScaleDown;
        a1( Rescale ) = a1( Rescale ) * ScaleDown;
        a2( Rescale ) = a2( Rescale ) * ScaleDown;
        a3( Rescale ) = a3( Rescale ) * ScaleDown;
    end
    q0 = a0 ./ b0;
    [ r0, r1, r2, r3 ] = DDTimesUnderlyingAsQD( b0, b1, q0 );
    [ r0, r1, r2, r3 ] = QDPlusQD( a0, a1, a2, a3, -r0, -r1, -r2, -r3 );
    q1 = r0 ./ b0;
    [ t0, t1, t2, t3 ] = DDTimesUnderlyingAsQD( b0, b1, q1 );
    [ r0, r1, r2, r3 ] = QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
    q2 = r0 ./ b0;
    [ t0, t1, t2, t3 ] = DDTimesUnderlyingAsQD( b0, b1, q2 );
    [ r0, r1, r2, r3 ] = QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
    q3 = r0 ./ b0;
    [ t0, t1, t2, t3 ] = DDTimesUnderlyingAsQD( b0, b1, q3 );
    [ r0, ~, ~, ~ ] = QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
    q4 = r0 ./ b0;
    [ s0, s1, s2, s3 ] = Renorm5( q0, q1, q2, q3, q4 );
    if any( Rescale, 'all' )
        ScaleUp = 2 ^ 53;
        s0( Rescale ) = s0( Rescale ) * ScaleUp;
        s1( Rescale ) = s1( Rescale ) * ScaleUp;
        s2( Rescale ) = s2( Rescale ) * ScaleUp;
        s3( Rescale ) = s3( Rescale ) * ScaleUp;
    end
end
