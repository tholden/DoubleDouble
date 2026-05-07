function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
    if ~isa( a, 'ED.BaseDoubleDouble' ) || ~isa( b, 'ED.BaseDoubleDouble' )
        % Standard Dekker TwoProd (works for doubles and mixed types)
        p1 = a .* b;
        [ a1, a2 ] = Split( a );
        [ b1, b2 ] = Split( b );
        [ t1, e1 ] = UnderlyingPlusUnderlying( a1 .* b1, -p1 );
        [ t2, e2 ] = UnderlyingPlusUnderlying( t1, a1 .* b2 );
        [ t3, e3 ] = UnderlyingPlusUnderlying( t2, a2 .* b1 );
        p2 = t3 + (a2 .* b2 + e1 + e2 + e3);
        return
    end

    % Both inputs are DoubleDoubles: exact decomposition to scalar doubles.
    % a = (a0 + a1), b = (b0 + b1) where a0,a1,b0,b1 are doubles.
    % Compute exact product as 4 TwoProds on doubles:
    %   a*b = a0*b0 + a0*b1 + a1*b0 + a1*b1
    % Each TwoProd gives (hi, lo), yielding 8 doubles total.
    % Accumulate into two DoubleDoubles [p1, p2].

    a0 = a.v1;  a1 = a.v2;
    b0 = b.v1;  b1 = b.v2;

    % 4 exact TwoProds on doubles (recursive call hits the Dekker path)
    [ h0, l0 ] = UnderlyingTimesUnderlying( a0, b0 );  % a0*b0
    [ h1, l1 ] = UnderlyingTimesUnderlying( a0, b1 );  % a0*b1
    [ h2, l2 ] = UnderlyingTimesUnderlying( a1, b0 );  % a1*b0
    [ h3, l3 ] = UnderlyingTimesUnderlying( a1, b1 );  % a1*b1

    % Accumulate 8 doubles into [d0,d1,d2,d3] (= p1 + p2 as two DDs)
    % using exact TwoSum cascades (same structure as DDTimesDDAsQD).

    % Tier 1: combine cross-term highs with leading error
    [ h1, h2, l0 ] = ThreeSum( h1, h2, l0 );

    % Tier 2: combine next-order terms
    [ h2, h3 ] = UnderlyingPlusUnderlying( h2, h3 );
    [ l1, l2 ] = UnderlyingPlusUnderlying( l1, l2 );

    [ s0, t0 ] = UnderlyingPlusUnderlying( h2, l1 );
    [ s1, t1 ] = UnderlyingPlusUnderlying( h3, l2 );
    [ s1, t0 ] = UnderlyingPlusUnderlying( s1, t0 );
    [ s2, e_s2 ] = UnderlyingPlusUnderlying( t0, t1 );

    p2val = s0;
    [ p3val, q0val, e_q0 ] = ThreeSum( l3, l0, s1 );
    [ p4val, e_p4 ] = UnderlyingPlusUnderlying( q0val, s2 );

    % Renorm5 into 4 non-overlapping doubles
    [ d0, d1, d2, d3 ] = Renorm5( h0, h1, p2val, p3val, p4val );

    % Add back the small error terms
    [ d0, d1, d2, d3 ] = QDPlusUnderlying( d0, d1, d2, d3, e_s2 );
    [ d0, d1, d2, d3 ] = QDPlusUnderlying( d0, d1, d2, d3, e_q0 );
    [ d0, d1, d2, d3 ] = QDPlusUnderlying( d0, d1, d2, d3, e_p4 );

    % Pack into two DoubleDoubles: p1 = (d0,d1), p2 = (d2,d3)
    p1 = a.Make( d0, d1 );
    p2 = a.Make( d2, d3 );
end
