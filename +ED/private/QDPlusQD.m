function [ x0, x1, x2, x3 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % IEEE-accurate QD+QD addition.
    % Uses pairwise TwoSum on matching components, then cascades
    % the error terms through ThreeSum to track all rounding errors.
    % Only the final lowest-order accumulation is sloppy.

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

    % Accumulate remaining residuals: e0, e1, e2, e3
    [ e0, e1_2 ] = UnderlyingPlusUnderlying( e0, e1 );
    [ e0, e1_3 ] = UnderlyingPlusUnderlying( e0, e2 );
    [ e0, e1_4 ] = UnderlyingPlusUnderlying( e0, e3 );
    e0 = e0 + ( e1_2 + e1_3 + e1_4 );

    [ x0, x1, x2, x3 ] = Renorm5( s0, s1, s2, s3, e0 );
end
