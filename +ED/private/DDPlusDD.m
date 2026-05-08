function [ s1, s2 ] = DDPlusDD( a1, a2, b1, b2 )
    % IEEE-accurate DD+DD addition ( cf. dd_real::ieee_add in QD library ).
    % Variable names match C++ where possible: s1,s2,t1,t2.
    % Improvement over C++: uses two_sum instead of sloppy s2+=t1.

    % if ~isreal( a1 ) || ~isreal( a2 ) || ~isreal( b1 ) || ~isreal( b2 )
    %     [ s1r, s2r ] = DDPlusDD( real( a1 ), real( a2 ), real( b1 ), real( b2 ) );
    %     [ s1i, s2i ] = DDPlusDD( imag( a1 ), imag( a2 ), imag( b1 ), imag( b2 ) );
    %     s1 = complex( s1r, s1i );
    %     s2 = complex( s2r, s2i );
    %     return
    % end

    % if abs( a2 ) > abs( b1 )
    %     a2_ = a2;
    %     a2 = b1;
    %     b1 = a2_;
    % end

    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ t1, t2 ] = UnderlyingPlusUnderlying( a2, b2 );
    [ s2, e ] = UnderlyingPlusUnderlying( s2, t1 );    % C++ does sloppy s2 += t1
    [ s1, s2 ] = DDNormalize( s1, s2 );
    tmp = t2 + e;
    s2 = s2 + tmp;
    [ s1, s2 ] = DDNormalize( s1, s2 );

end
