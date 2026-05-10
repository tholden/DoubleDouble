function [ x0, x1, x2, x3, x4, x5, x6, x7 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
    % IEEE-accurate QD+QD addition.
    % Sorts all 8 components by decreasing magnitude using a manual Sort8 network, 
    % then cascades via exact two_sum accumulation returning all 8 exact components.

    if ~isreal( a0 ) || ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( a3 ) || ...
            ~isreal( b0 ) || ~isreal( b1 ) || ~isreal( b2 ) || ~isreal( b3 )
        [ x0r, x1r, x2r, x3r, x4r, x5r, x6r, x7r ] = QDPlusQD( real( a0 ), real( a1 ), real( a2 ), real( a3 ),  ...
            real( b0 ), real( b1 ), real( b2 ), real( b3 ) );
        [ x0i, x1i, x2i, x3i, x4i, x5i, x6i, x7i ] = QDPlusQD( imag( a0 ), imag( a1 ), imag( a2 ), imag( a3 ),  ...
            imag( b0 ), imag( b1 ), imag( b2 ), imag( b3 ) );
        x0 = complex( x0r, x0i );
        x1 = complex( x1r, x1i );
        x2 = complex( x2r, x2i );
        x3 = complex( x3r, x3i );
        x4 = complex( x4r, x4i );
        x5 = complex( x5r, x5i );
        x6 = complex( x6r, x6i );
        x7 = complex( x7r, x7i );
        return
    end

    [ a0, a1, a2, a3, b0, b1, b2, b3 ] = ED.ExpandSingleton( a0, a1, a2, a3, b0, b1, b2, b3 );
    
    [ z1, z2, z3, z4, z5, z6, z7, z8 ] = Sort8( a0, a1, a2, a3, b0, b1, b2, b3 );

    [ s1, s2 ] = UnderlyingPlusUnderlying( z1, z2 );
    [ s2, s3 ] = UnderlyingPlusUnderlying( s2, z3 );
    [ s3, s4 ] = UnderlyingPlusUnderlying( s3, z4 );
    [ s4, s5 ] = UnderlyingPlusUnderlying( s4, z5 );
    [ s5, s6 ] = UnderlyingPlusUnderlying( s5, z6 );
    [ s6, s7 ] = UnderlyingPlusUnderlying( s6, z7 );
    [ s7, s8 ] = UnderlyingPlusUnderlying( s7, z8 );

    [ s1, s2, s3, s4 ] = QDNormalize( s1, s2, s3, s4 );
    
    % Proper renormalization of the bottom 4 components
    % We must cascade upwards from s8
    [ s5, s6, s7, s8 ] = QDNormalize( s5, s6, s7, s8 );

    % Now we have a normalized [s1..s4] block and a normalized [s5..s8] block.
    % Wait! We don't want them separated, they might overlap!
    % If we just return them, BaseQuadDouble will package them.
    x0 = s1; x1 = s2; x2 = s3; x3 = s4;
    x4 = s5; x5 = s6; x6 = s7; x7 = s8;
    
end

function [ a, b, c, d, e, f, g, h ] = Sort8( a, b, c, d, e, f, g, h )
    aa = abs( a ); ab = abs( b ); ac = abs( c ); ad = abs( d );
    ae = abs( e ); af = abs( f ); ag = abs( g ); ah = abs( h );
    
    [a, b, aa, ab] = swap_if(a, b, aa, ab); [c, d, ac, ad] = swap_if(c, d, ac, ad); [e, f, ae, af] = swap_if(e, f, ae, af); [g, h, ag, ah] = swap_if(g, h, ag, ah); 
    [a, c, aa, ac] = swap_if(a, c, aa, ac); [b, d, ab, ad] = swap_if(b, d, ab, ad); [e, g, ae, ag] = swap_if(e, g, ae, ag); [f, h, af, ah] = swap_if(f, h, af, ah); 
    [b, c, ab, ac] = swap_if(b, c, ab, ac); [f, g, af, ag] = swap_if(f, g, af, ag); [a, e, aa, ae] = swap_if(a, e, aa, ae); [b, f, ab, af] = swap_if(b, f, ab, af); 
    [c, g, ac, ag] = swap_if(c, g, ac, ag); [d, h, ad, ah] = swap_if(d, h, ad, ah); [c, e, ac, ae] = swap_if(c, e, ac, ae); [d, f, ad, af] = swap_if(d, f, ad, af); 
    [b, c, ab, ac] = swap_if(b, c, ab, ac); [d, e, ad, ae] = swap_if(d, e, ad, ae); [f, g, af, ag] = swap_if(f, g, af, ag); 
end

function [ v1, v2, a1, a2 ] = swap_if( v1, v2, a1, a2 )
    swap = a1 < a2;
    if any( swap, 'all' )
        tmp = v1(swap); v1(swap) = v2(swap); v2(swap) = tmp;
        tmp = a1(swap); a1(swap) = a2(swap); a2(swap) = tmp;
    end
end
