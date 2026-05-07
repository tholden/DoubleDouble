function [ varargout ] = ExpandSingleton( varargin )
    l = cellfun( @( x ) length( size( x ) ), varargin, 'UniformOutput', true );
    n = length( varargin );
    ss = ones( n, max( l ) );
    for i = 1 : n
        ss( i, 1 : l( i ) ) = size( varargin{ i } );
    end
    if all( ss == ss( 1, : ), 'all' )
        varargout = varargin;
        return
    end
    s = max( ss, [], 1 );
    isZeroDim = any( ss == 0, 1 );
    if any( isZeroDim )
        if any( any( ss( :, isZeroDim ) > 1 ) )
            error( 'Arrays have incompatible sizes for this operation.' );
        end
        s( isZeroDim ) = 0;
    end
    if all( s <= 1 )
        varargout = varargin;
        return
    end
    varargout = cell( 1, n );
    for i = 1 : n
        RepFactor = ones( 1, max( length( s ), 2 ) );
        mask = find( ( ss( i, : ) == 1 ) & ( s ~= 1 ) );
        RepFactor( mask ) = s( mask );
        varargout{ i } = repmat( varargin{ i }, RepFactor );
    end
end
