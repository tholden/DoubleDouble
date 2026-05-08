function [ a1, a2 ] = Split( a )

    if isreal( a )
        Select = ( a > 6.6969287949141707559e+299 ) | ( a < -6.6969287949141707559e+299 ); % 2^996
        a( Select ) = a( Select ) * 3.7252902984619140625e-09; % 2^( -28 )
        t1 = 134217729.0 * a; % 2^27 + 1
        t2 = t1 - a;
        a1 = t1 - t2;
        a2 = a - a1;
        a1( Select ) = a1( Select ) * 268435456.0; % 2^28
        a2( Select ) = a2( Select ) * 268435456.0; % 2^28
    else
        [ r1, r2 ] = Split( real( a ) );
        [ i1, i2 ] = Split( imag( a ) );
        a1 = complex( r1, i1 );
        a2 = complex( r2, i2 );
    end

end
