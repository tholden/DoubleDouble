function [ s1, s2, s3, s4 ] = DDPlusDD( a1, a2, b1, b2 )
    % IEEE-accurate DD+DD addition.
    % Utilizes a manual sorting network to order the 4 components by absolute magnitude.
    % Returns 4 exact components.

    if ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( b1 ) || ~isreal( b2 )
        [ s1r, s2r, s3r, s4r ] = DDPlusDD( real( a1 ), real( a2 ), real( b1 ), real( b2 ) );
        [ s1i, s2i, s3i, s4i ] = DDPlusDD( imag( a1 ), imag( a2 ), imag( b1 ), imag( b2 ) );
        s1 = complex( s1r, s1i );
        s2 = complex( s2r, s2i );
        s3 = complex( s3r, s3i );
        s4 = complex( s4r, s4i );
        return
    end

    [ a1, a2, b1, b2 ] = ED.ExpandSingleton( a1, a2, b1, b2 );

    [ z1, z2, z3, z4 ] = Sort4( a1, a2, b1, b2 );

    [ s1, s2 ] = UnderlyingPlusUnderlying( z1, z2 );
    [ s2, s3 ] = UnderlyingPlusUnderlying( s2, z3 );
    [ s3, s4 ] = UnderlyingPlusUnderlying( s3, z4 );
    
    [ s1, s2, s3, s4 ] = QDNormalize( s1, s2, s3, s4 );

end

function [ a, b, c, d ] = Sort4( a, b, c, d )
    aa = abs( a ); ab = abs( b ); ac = abs( c ); ad = abs( d );
    
    % Layer 1
    swap = aa < ab;
    if any( swap, 'all' )
        tmp = a(swap); a(swap) = b(swap); b(swap) = tmp;
        tmp = aa(swap); aa(swap) = ab(swap); ab(swap) = tmp;
    end
    
    swap = ac < ad;
    if any( swap, 'all' )
        tmp = c(swap); c(swap) = d(swap); d(swap) = tmp;
        tmp = ac(swap); ac(swap) = ad(swap); ad(swap) = tmp;
    end
    
    % Layer 2
    swap = aa < ac;
    if any( swap, 'all' )
        tmp = a(swap); a(swap) = c(swap); c(swap) = tmp;
        tmp = aa(swap); aa(swap) = ac(swap); ac(swap) = tmp;
    end
    
    swap = ab < ad;
    if any( swap, 'all' )
        tmp = b(swap); b(swap) = d(swap); d(swap) = tmp;
        tmp = ab(swap); ab(swap) = ad(swap); ad(swap) = tmp;
    end
    
    % Layer 3
    swap = ab < ac;
    if any( swap, 'all' )
        tmp = b(swap); b(swap) = c(swap); c(swap) = tmp;
        % tmp = ab(swap); ab(swap) = ac(swap); ac(swap) = tmp;
    end
end
