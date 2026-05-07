function [ x0, x1, x2, x3 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % IEEE-accurate QD+QD addition ( cf. qd_real::ieee_add in QD library ).
    % Sorts all 8 components by decreasing magnitude, then cascades via
    % the sloppy_add accumulation structure for vectorized operation.

    [ a0, a1, a2, a3, b0, b1, b2, b3 ] = ED.ExpandSingleton( a0, a1, a2, a3, b0, b1, b2, b3 );
    CatDim = ndims( a0 ) + 1;
    z = cat( CatDim, a0, b0, a1, b1, a2, b2, a3, b3 );

    if isreal( z )
        z = sort( z, CatDim, 'descend', 'ComparisonMethod', 'abs' );
    else
        zr = real( z );
        zi = imag( z );
        zr = sort( zr, CatDim, 'descend', 'ComparisonMethod', 'abs' );
        zi = sort( zi, CatDim, 'descend', 'ComparisonMethod', 'abs' );
        z = complex( zr, zi );
    end

    sz = size( z );
    zFlat = reshape( z, [], sz( end ) );
    a0 = reshape( zFlat( :, 1 ), sz( 1 : ( end - 1 ) ) );
    b0 = reshape( zFlat( :, 2 ), sz( 1 : ( end - 1 ) ) );
    a1 = reshape( zFlat( :, 3 ), sz( 1 : ( end - 1 ) ) );
    b1 = reshape( zFlat( :, 4 ), sz( 1 : ( end - 1 ) ) );
    a2 = reshape( zFlat( :, 5 ), sz( 1 : ( end - 1 ) ) );
    b2 = reshape( zFlat( :, 6 ), sz( 1 : ( end - 1 ) ) );
    a3 = reshape( zFlat( :, 7 ), sz( 1 : ( end - 1 ) ) );
    b3 = reshape( zFlat( :, 8 ), sz( 1 : ( end - 1 ) ) );

    % Sloppy-add cascade on the sorted components.
    [ s0, t0 ] = UnderlyingPlusUnderlying( a0, b0 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ s2, t2 ] = UnderlyingPlusUnderlying( a2, b2 );
    [ s3, t3 ] = UnderlyingPlusUnderlying( a3, b3 );

    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    [ s2, t0, t1 ] = ThreeSum( s2, t0, t1 );
    [ s3, t0 ] = ThreeSum2( s3, t0, t2 );
    t0 = t0 + t1 + t3;

    [ x0, x1, x2, x3 ] = Renorm5( s0, s1, s2, s3, t0 );
end
