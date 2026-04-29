function [ v, d ] = eigExt( x )

    [ v, d ] = eig( x.v1 );
    v = x.Promote( v );
    d = x.Promote( diag( d ) );
    C = length( d );
    I = eye( C, 'like', x );

    for c = 1 : C

        vi = v( :, c );
        dii = d( c, 1 );
        err = Inf;

        while true
            nvi = ( x - dii * I ) \ vi;
            nvi = nvi ./ norm( nvi );
            if any( ~isfinite( nvi ) )
                break
            end
            vi = nvi;
            odii = dii;
            xTvi = x * vi;
            dii = ( vi' * xTvi ) ./ ( vi' * vi );
            oerr = err;
            errv = abs( xTvi - dii * vi );
            err = sum( errv .* errv );
            if err > oerr
                dii = odii;
                break
            end
            if ( err == 0 ) || ( err == oerr )
                break
            end
        end

        d( c, 1 ) = dii;
        v( :, c ) = vi;

    end

    if nargout < 2
        v = d;
    else
        d = diag( d );
    end

end
