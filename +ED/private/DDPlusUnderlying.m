function [ s1, s2 ] = DDPlusUnderlying( a1, a2, b )
    % DD + underlying addition ( cf. operator+(dd_real, double) in QD library ).
    % Variable names match C++: s1,s2.
    % Improvement over C++: uses two_sum instead of sloppy s2 += a[1].
    [ s1, s2 ] = UnderlyingPlusUnderlying( a1, b );
    [ s2, e ] = UnderlyingPlusUnderlying( s2, a2 );
    [ s1, s2 ] = DDNormalize( s1, s2 );
    s2 = s2 + e;
    % Paranoid alternative to s2 = s2 + e:
    % [ s2, e2 ] = UnderlyingPlusUnderlying( s2, e );
    % [ s1, s2 ] = DDNormalize( s1, s2 );
    % s2 = s2 + e2;
    [ s1, s2 ] = DDNormalize( s1, s2 );
end
