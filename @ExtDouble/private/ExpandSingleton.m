function [ varargout ] = ExpandSingleton( varargin )
    l = cellfun( @( x ) length( size( x ) ), varargin, 'UniformOutput', true );
    n = length( varargin );
    ss = ones( n, max( l ) );
    for i = 1 : n
        ss( i, 1 : l( i ) ) = size( varargin{ i } );
    end
    s = max( ss, [], 1 );
    isZeroDim = any( ss == 0, 1 );
    if any( isZeroDim )
        if any( any( ss( :, isZeroDim ) > 1 ) )
            error( 'Arrays have incompatible sizes for this operation.' );
        end
        s( isZeroDim ) = 0;
    end
    varargout = cell( 1, n );
    for i = 1 : n
        repfac = ones( 1, length( s ) );
        mask = ( ss( i, : ) == 1 ) & ( s ~= 1 );
        repfac( mask ) = s( mask );
        varargout{ i } = repmat( varargin{ i }, repfac );
    end
end
