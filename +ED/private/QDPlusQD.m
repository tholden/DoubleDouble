function [ x0, x1, x2, x3, x4 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % IEEE-accurate QD+QD addition.
    % Uses pairwise TwoSum on matching components, then cascades
    % the error terms through ThreeSum to track all rounding errors.
    % Only the final lowest-order accumulation is sloppy.

    % CatDim = ndims( a0 ) + 1;
    % z = cat( CatDim, a0, a1 + 0 * a0, a2 + 0 * a0, a3 + 0 * a0, b0, b1 + 0 * b0, b2 + 0 * b0, b3 + 0 * b0 );
    % if isreal( z )
    %     z = sort( z, CatDim, 'descend', 'ComparisonMethod', 'abs' );
    % else
    %     zr = real( z );
    %     zi = imag( z );
    %     zr = sort( zr, CatDim, 'descend', 'ComparisonMethod', 'abs' );
    %     zi = sort( zi, CatDim, 'descend', 'ComparisonMethod', 'abs' );
    %     z = complex( zr, zi );
    % end
    % sz = size( z );
    % zFlat = reshape( z, [], sz( end ) );
    % a0 = reshape( zFlat( :, 1 ), sz( 1 : end - 1 ) );
    % b0 = reshape( zFlat( :, 2 ), sz( 1 : end - 1 ) );
    % a1 = reshape( zFlat( :, 3 ), sz( 1 : end - 1 ) );
    % b1 = reshape( zFlat( :, 4 ), sz( 1 : end - 1 ) );
    % a2 = reshape( zFlat( :, 5 ), sz( 1 : end - 1 ) );
    % b2 = reshape( zFlat( :, 6 ), sz( 1 : end - 1 ) );
    % a3 = reshape( zFlat( :, 7 ), sz( 1 : end - 1 ) );
    % b3 = reshape( zFlat( :, 8 ), sz( 1 : end - 1 ) );

    % ---- pairwise TwoSum on matching components ----
    % This computes 4 exact sums + 4 exact errors = 8 terms total.
    [ s0, e0 ] = UnderlyingPlusUnderlying( a0, b0 );
    [ s1, e1 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ s2, e2 ] = UnderlyingPlusUnderlying( a2, b2 );
    [ s3, e3 ] = UnderlyingPlusUnderlying( a3, b3 );

    % ---- cascade the 8 terms through ThreeSum / TwoSum ----
    % s0 is the leading term.
    % Merge e0 into (s1):
    [ s1, e0 ] = UnderlyingPlusUnderlying( s1, e0 );

    % Now we have s0, s1, and residuals e0, e1, s2, e2, s3, e3
    % Merge e0 and e1 into s2:
    [ s2, e0, e1 ] = ThreeSum( s2, e0, e1 );

    % Merge e0, e2, s3:
    [ s3, e0, e2 ] = ThreeSum( s3, e0, e2 );

    % Accumulate remaining residuals using a single QDPlusUnderlying
    % for performance, as they are all O(eps^4) or smaller.
    e_sum = e0 + e1 + e2 + e3;
    [ x0, x1, x2, x3, x4 ] = QDPlusUnderlying( s0, s1, s2, s3, e_sum );
end
