function [ s1, s2, s3, s4, s5, s6, s7, s8 ] = QDPlusUnderlying( a1, a2, a3, a4, b )
    % QD + underlying addition.
    % Returns 8 exact components.

    if ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( a3 ) || ~isreal( a4 ) || ~isreal( b )
        [ s1r, s2r, s3r, s4r, s5r, s6r, s7r, s8r ] = QDPlusUnderlying( real( a1 ), real( a2 ), real( a3 ), real( a4 ), real( b ) );
        [ s1i, s2i, s3i, s4i, s5i, s6i, s7i, s8i ] = QDPlusUnderlying( imag( a1 ), imag( a2 ), imag( a3 ), imag( a4 ), imag( b ) );
        s1 = complex( s1r, s1i );
        s2 = complex( s2r, s2i );
        s3 = complex( s3r, s3i );
        s4 = complex( s4r, s4i );
        s5 = complex( s5r, s5i );
        s6 = complex( s6r, s6i );
        s7 = complex( s7r, s7i );
        s8 = complex( s8r, s8i );
        return
    end

    b2 = zeros( size( b ), 'like', b );
    b3 = zeros( size( b ), 'like', b );
    b4 = zeros( size( b ), 'like', b );
    
    [ s1, s2, s3, s4, s5, s6, s7, s8 ] = QDPlusQD( a1, a2, a3, a4, b, b2, b3, b4 );
end
