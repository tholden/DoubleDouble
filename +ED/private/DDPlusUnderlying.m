function [ s1, s2, s3, s4 ] = DDPlusUnderlying( a1, a2, b )
    % DD + underlying addition.
    % Returns 4 exact components.

    if ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( b )
        [ s1r, s2r, s3r, s4r ] = DDPlusUnderlying( real( a1 ), real( a2 ), real( b ) );
        [ s1i, s2i, s3i, s4i ] = DDPlusUnderlying( imag( a1 ), imag( a2 ), imag( b ) );
        s1 = complex( s1r, s1i );
        s2 = complex( s2r, s2i );
        s3 = complex( s3r, s3i );
        s4 = complex( s4r, s4i );
        return
    end

    b2 = zeros( size( b ), 'like', b );
    [ s1, s2, s3, s4 ] = DDPlusDD( a1, a2, b, b2 );

end
