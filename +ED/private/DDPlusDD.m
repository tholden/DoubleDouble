function [ x1, x2 ] = DDPlusDD( a1, a2, b1, b2 )
    % IEEE-accurate DD+DD addition.
    % Sorts all 4 components by decreasing magnitude, then cascades via
    % two_sum accumulation (matching the QDPlusQD sloppy_add approach).

    if ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( b1 ) || ~isreal( b2 )
        [ x1r, x2r ] = DDPlusDD( real( a1 ), real( a2 ), real( b1 ), real( b2 ) );
        [ x1i, x2i ] = DDPlusDD( imag( a1 ), imag( a2 ), imag( b1 ), imag( b2 ) );
        x1 = complex( x1r, x1i );
        x2 = complex( x2r, x2i );
        return
    end

    [ a1, a2, b1, b2 ] = ED.ExpandSingleton( a1, a2, b1, b2 );

    % Sort 4 components by |value| descending using a sorting network.
    % For compound types, use double(abs()) for comparisons.
    c1 = a1; c2 = b1; c3 = a2; c4 = b2;

    % 4-element sorting network (5 compare-and-swap operations)
    % Pairs: (1,2), (3,4), (1,3), (2,4), (2,3)
    [ c1, c2 ] = CompareSwap( c1, c2 );
    [ c3, c4 ] = CompareSwap( c3, c4 );
    [ c1, c3 ] = CompareSwap( c1, c3 );
    [ c2, c4 ] = CompareSwap( c2, c4 );
    [ c2, c3 ] = CompareSwap( c2, c3 );

    u = c1;
    v = c2;

    % two_sum(u, v)
    s_uv = u + v;
    bb = s_uv - u;
    tmp1 = s_uv - bb;
    tmp2 = u - tmp1;
    tmp3 = v - bb;
    v = tmp2 + tmp3;
    u = s_uv;

    % Process c3
    t = c3;
    [ u, v ] = AccumulateOne( u, v, t );

    % Process c4
    t = c4;
    [ u, v ] = AccumulateOne( u, v, t );

    [ x1, x2 ] = DDNormalize( u, v );
end

function [ u, v ] = AccumulateOne( u, v, t )
    % two_sum(v, t)
    s1 = v + t;
    bb1 = s1 - v;
    tmp1 = s1 - bb1;
    tmp2 = v - tmp1;
    tmp3 = t - bb1;
    v_err = tmp2 + tmp3;

    % two_sum(u, s1)
    s2 = u + s1;
    bb2 = s2 - u;
    tmp1 = s2 - bb2;
    tmp2 = u - tmp1;
    tmp3 = s1 - bb2;
    u_err = tmp2 + tmp3;

    % Fold corrections back into v
    v = v_err + u_err;
    u = s2;
end

function [ a, b ] = CompareSwap( a, b )
    % Sort a, b so |a| >= |b|, element-wise
    swap = abs( a ) < abs( b );
    if any( swap, 'all' )
        a_old = a;
        a( swap ) = b( swap );
        b( swap ) = a_old( swap );
    end
end
