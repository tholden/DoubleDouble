function [ s1, s2 ] = DDPlusDD( a1, a2, b1, b2 )
    % IEEE-accurate DD+DD addition ( cf. dd_real::ieee_add in QD library ).
    % Variable names match C++ where possible: s1,s2,t1,t2.
    % Improvement over C++: uses two_sum instead of sloppy s2+=t1.
    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b1 );
    [ t1, t2 ] = UnderlyingPlusUnderlying( a2, b2 );
    [ s2, e ] = UnderlyingPlusUnderlying( s2, t1 );    % C++ does sloppy s2 += t1
    [ s1, s2 ] = DDNormalize( s1, s2 );
    s2 = s2 + (t2 + e);
    [ s1, s2 ] = DDNormalize( s1, s2 );
end
