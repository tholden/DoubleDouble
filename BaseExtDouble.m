classdef (Abstract) BaseExtDouble

    properties ( SetAccess = public, GetAccess = public )
        v1
        v2
    end

    properties ( Constant, GetAccess = public )
        SingletonExpansionNotSupported = ~BaseExtDouble.TestSingletonExpansion();
    end

    properties ( Abstract, Constant )
        zero
        one
        tiny
        pi
        InverseFactorial
        NInverseFactorial
        piT2
        piD2
        piD16
        log_2
        log_10
        SinTable
        CosTable
        ExpRescale
        LogSteps
    end

    methods (Abstract)
        v = Make( z, a1, a2 )
        v = Promote( v, a )
    end

    methods (Sealed)

        function v = cat( dim, a, varargin )
            if isempty( varargin )
                v = a;
            else
                b = varargin{ 1 };
                if length( varargin ) > 1
                    v = cat( dim, cat( dim, a, b ), varargin{ 2 : end } );
                else
                    if ~isa( a, 'BaseExtDouble' )
                        a = b.Promote( a );
                    end
                    if ~isa( b, 'BaseExtDouble' )
                        b = a.Promote( b );
                    end
                    x1 = cat( dim, a.v1, b.v1 );
                    x2 = cat( dim, a.v2, b.v2 );
                    v = a.Make( x1, x2 );
                end
            end
        end

        function v = horzcat( a, varargin )
            if isempty( varargin )
                v = a;
            else
                b = varargin{ 1 };
                if length( varargin ) > 1
                    v = cat( 2, cat( 2, a, b ), varargin{ 2 : end } );
                else
                    v = cat( 2, a, b );
                end
            end
        end

        function v = vertcat( a, varargin )
            if isempty( varargin )
                v = a;
            else
                b = varargin{ 1 };
                if length( varargin ) > 1
                    v = cat( 1, cat( 1, a, b ), varargin{ 2 : end } );
                else
                    v = cat( 1, a, b );
                end
            end
        end

        function v = Index( v, varargin )
            assert( isa( v, 'BaseExtDouble' ) );
            v.v1 = Index( v.v1, varargin{:} );
            v.v2 = Index( v.v2, varargin{:} );
        end

        function v = Assign( v, Value, varargin )
            if ~isa( v, 'BaseExtDouble' )
                assert( isa( Value, 'BaseExtDouble' ) );
                if isempty( v )
                    v = zeros( cellfun( @max, varargin, 'UniformOutput', true ), 'like', Value );
                else
                    v = Value.Promote( v );
                end
            elseif ~isa( Value, 'BaseExtDouble' )
                Value = v.Promote( Value );
            end
            v.v1 = Assign( v.v1, Value.v1, varargin{:} );
            v.v2 = Assign( v.v2, Value.v2, varargin{:} );
        end



        function varargout = ToSumOfDoubles( v )
            varargout = cell( 1, nargout );
            [ varargout{:} ] = ToSumOfDoubles( v.v1 );
            FirstEmpty = find( cellfun( @isempty, varargout, 'UniformOutput', true ), 1, 'first' );
            if ~isempty( FirstEmpty )
                [ varargout{ FirstEmpty : nargout } ] = ToSumOfDoubles( v.v2 );
            end
        end


        function disp( v )
            if isempty( v.v1 )
                disp( '     []' );
            else
                disp( v.v1 );
            end
            disp( '     +' );
            if isempty( v.v2 )
                disp( '     []' );
            else
                disp( v.v2 );
            end
            disp( ' ' );
        end

        function v = double( v )
            v = double( v.v1 );
        end

        function v = isscalar( v )
            v = isscalar( v.v1 );
        end

        function v = isreal( v )
            v = isreal( v.v1 ) && isreal( v.v2 );
        end

        function v = isnumeric( v )
            v = isnumeric( v.v1 ) && isnumeric( v.v2 );
        end

        function v = isfinite( v )
            v = isfinite( v.v1 ) & isfinite( v.v2 );
        end

        function v = isinf( v )
            v = isinf( v.v1 ) | isinf( v.v2 );
        end

        function v = isnan( v )
            v = isnan( v.v1 ) | isnan( v.v2 );
        end

        function v = sparse( i, j, v, m, n, nz )
            if nargin == 0
                v = i.Make( sparse( [] ), sparse( [] ) );
            elseif nargin < 3
                assert( nargin == 1 );
                v = i;
                v.v1 = sparse( v.v1 );
                v.v2 = sparse( v.v2 );
            elseif nargin == 3
                v.v1 = sparse( i, j, v.v1 );
                v.v2 = sparse( i, j, v.v2 );
            elseif nargin == 4
                v.v1 = sparse( i, j, v.v1, m, n );
                v.v2 = sparse( i, j, v.v2, m, n );
            else
                v.v1 = sparse( i, j, v.v1, m, n, nz );
                v.v2 = sparse( i, j, v.v2, m, n, nz );
            end
        end

        function v = any( v, varargin )
            v = any( v ~= 0, varargin{:} );
        end

        function v = all( v, varargin )
            v = all( v ~= 0, varargin{:} );
        end

        function varargout = find( v, varargin )
            if nargout == 1
                varargout{ 1 } = find( v ~= 0, varargin{:} );
            elseif nargout >= 2
                [ varargout{ 1 }, varargout{ 2 } ] = find( v ~= 0, varargin{:} );
                if nargout >= 3
                    LinearIndex = sub2ind( size( v ), varargout{ 1 }, varargout{ 2 } );
                    varargout{ 3 } = v.Make( Index( v.v1, LinearIndex ), Index( v.v2, LinearIndex ) );
                end
            end
        end

        function v = real( v )
            v.v1 = real( v.v1 );
            v.v2 = real( v.v2 );
        end

        function v = imag( v )
            v.v1 = imag( v.v1 );
            v.v2 = imag( v.v2 );
        end

        function v = conj( v )
            v.v1 = conj( v.v1 );
            v.v2 = conj( v.v2 );
        end

        function v = angle( v )
            if isreal( v )
                Select = v >= 0;
                v.v1 = Assign( v.v1, 0, Select );
                v.v2 = Assign( v.v2, 0, Select );
                v.v1 = Assign( v.v1, v.pi.v1, ~Select );
                v.v2 = Assign( v.v2, v.pi.v2, ~Select );
            else
                v = atan2( imag( v ), real( v ) );
            end
        end

        function [ v, varargout ] = size( v, varargin )
            v = size( v.v1, varargin{:} );
            if nargout > 1
                varargout = num2cell( v( 2:end ) );
                v = v( 1 );
            end
        end

        function v = length( v )
            v = max( size( v ) );
        end

        function v = numel( v )
            v = numel( v.v1 );
        end

        function n = numArgumentsFromSubscript( v, s, IndexingContext )
            if strcmp( s(1).type, '.' )
                n = builtin( 'numArgumentsFromSubscript', v, s, IndexingContext );
            else
                n = numArgumentsFromSubscript( v.v1, s, IndexingContext );
            end
        end

        function v = end( v, k, n )
            if n == 1
                v = numel( v.v1 );
            else
                v = size( v.v1 );
                if k <= length( v )
                    v = v( k );
                else
                    v = 1;
                end
            end
        end

        function v = repmat( v, varargin )
            v = v.Make( repmat( v.v1, varargin{:} ), repmat( v.v2, varargin{:} ) );
        end

        function v = reshape( v, varargin )
            v.v1 = reshape( v.v1, varargin{:} );
            v.v2 = reshape( v.v2, varargin{:} );
        end

        function v = isequal( a, b, varargin )
            if any( size( a ) ~= size( b ) )
                v = false;
                return
            end
            v = a == b;
            v = all( v(:) );
            if nargin > 2
                for i = 1 : length( varargin )
                    if ~v
                        break
                    end
                    v = v && isequal( a, varargin{i} );
                end
            end
        end

        function v = isempty( v )
            if builtin('numel', v) ~= 1
                v = builtin('isempty', v);
            else
                v = isempty( v.v1 );
            end
        end

        function v = diag( v, k )
            if nargin < 2
                v = v.Make( diag( v.v1 ), diag( v.v2 ) );
            else
                v = v.Make( diag( v.v1, k ), diag( v.v2, k ) );
            end
        end

        function v = tril( v, k )
            if nargin < 2
                v = v.Make( tril( v.v1 ), tril( v.v2 ) );
            else
                v = v.Make( tril( v.v1, k ), tril( v.v2, k ) );
            end
        end

        function v = triu( v, k )
            if nargin < 2
                v = v.Make( triu( v.v1 ), triu( v.v2 ) );
            else
                v = v.Make( triu( v.v1, k ), triu( v.v2, k ) );
            end
        end

        function v = plus( a, b )
            if isa( a, 'BaseExtDouble' )
                v = a.Plus( b );
            else
                v = b.Plus( a );
            end
        end

        function v = minus( a, b )
            if isa( a, 'BaseExtDouble' )
                v = a.Minus( b );
            else
                v = -( b.Minus( a ) );
            end
        end

        function v = uminus( v )
            v.v1 = -v.v1;
            v.v2 = -v.v2;
        end

        function v = uplus( v )
        end

        function v = times( a, b )
            if isa( a, 'BaseExtDouble' )
                v = a.Times( b );
            else
                v = b.Times( a );
            end
        end

        function v = mtimes( a, b )
            if ~isa( a, 'BaseExtDouble' )
                a = b.Promote( a );
            end
            if ~isa( b, 'BaseExtDouble' )
                b = a.Promote( b );
            end
            v = a.MTimes( b );
        end

        function v = rdivide( a, b )
            if isa( a, 'BaseExtDouble' )
                v = a.RDivide( b );
            else
                a = b.Promote( a );
                v = a.RDivide( b );
            end
        end

        function v = ldivide( a, b )
            if isa( b, 'BaseExtDouble' )
                v = b.RDivide( a );
            else
                b = a.Promote( b );
                v = b.RDivide( a );
            end
        end

        function v = mldivide( a, v )
            [ Ra, Ca ] = size( a );
            [ Rv, Cv ] = size( v );
            if ( ( Ra == 1 ) && ( Ca == 1 ) ) || ( ( Rv == 1 ) && ( Cv == 1 ) )
                v = a.LDivide( v );
                return
            end
            if ~isa( v, 'BaseExtDouble' )
                v = a.Promote( v );
            end
            assert( Ra == Rv );
            if Ra ~= Ca
                [ q, a ] = qr( a );
                v = q' * v;
                v = BackSubstitutionExt( v, a );
                return
            end
            if BaseExtDouble.IsEqualWithExpansion( triu( a, 1 ), 0 )
                % Lower triangular
                v = ForwardEliminationExt( v, a );
                return
            elseif BaseExtDouble.IsEqualWithExpansion( tril( a, -1 ), 0 )
                % Upper triangular
                v = BackSubstitutionExt( v, a );
                return
            elseif BaseExtDouble.IsEqualWithExpansion( a, a' )
                [ L, d ] = ldl( a, 'vector_d' );
                if all( all( isfinite( L ) ) ) && all( isfinite( d ) )
                    % Positive definite
                    v = ForwardEliminationExt( v, L );
                    v = v ./ d;
                    v = BackSubstitutionExt( v, L' );
                    return
                end
            end
            % Triangular factorization
            [ L, U, p ] = lu( a, 'vector' );

            % Permutation and forward elimination
            v.v1 = Index( v.v1, p, ':' );
            v.v2 = Index( v.v2, p, ':' );
            v = ForwardEliminationExt( v, L );

            % Back substitution
            v = BackSubstitutionExt( v, U );
        end

        function v = mrdivide( v, a )
            if ~isa( a, 'BaseExtDouble' )
                a = v.Promote( a );
            end
            if ~isa( v, 'BaseExtDouble' )
                v = a.Promote( v );
            end
            at = a';
            vt = v';
            v = MLDivide(at, vt )';
        end

        function v = power( a, b )
            if ~isa( a, 'BaseExtDouble' )
                a = b.Promote( a );
            end
            if isa( b, 'BaseExtDouble' )
                b2 = b.TimesPowerOf2( 2 );
            else
                b2 = 2 * b;
            end
            bHalfInteger = ( double( b2 ) == b2 ) & ( floor( b2 ) == b2 );
            if ~any( bHalfInteger )
                v = exp( b .* log( a ) );
            else
                [ a, b, bHalfInteger ] = BaseExtDouble.ExpandSingleton( a, b, bHalfInteger );
                v = a.ones( size( a ) );
                bNonHalfInteger = find( ~bHalfInteger );
                bHalfInteger = find( bHalfInteger );
                if ~isempty( bNonHalfInteger )
                    v =  v.Assign(exp(  b(bNonHalfInteger ) .* log(  a.Index(bNonHalfInteger ) ) ), bNonHalfInteger  );
                end
                a =  a.Index(bHalfInteger );
                b = double(  b(bHalfInteger ) );
                Select = find( b < 0 );
                a =  a.Assign(1 ./  a.Index(Select ), Select );
                b(Select) = -b( Select );
                vv =   v.Index(bHalfInteger );
                Select = find( b ~= floor( b ) );
                vv =  vv.Assign( sqrt(  a.Index(Select ) ), Select  );
                Binary = dec2bin( floor( b ) );
                N = size( Binary, 2 );
                Power = a;
                Select = find( Binary( :, end ) == '1' );
                vv =  vv.Assign( a.Index(Select ), Select );
                for n = 2 : N
                    Power = Power .* Power;
                    Select = find( Binary( :, end + 1 - n ) == '1' );
                    vv =  vv.Assign( vv.Index(Select ) .*  Power.Index(Select ), Select );
                end
                v =  v.Assign(vv, bHalfInteger );
            end
        end

        function v = mpower( a, b )
            na = numel( a );
            nb = numel( b );
            if na <= 1
                if nb <= 1
                    v = a .^ b;
                else
                    [ v, d ] = eig( b );
                    d = diag( d );
                    assert( length( unique( d.v1 ) ) == length( d.v1 ) );
                    v = v * diag( a .^ d ) / v;
                end
            else
                if nb == 1
                    [ v, d ] = eig( a );
                    d = diag( d );
                    assert( length( unique( d.v1 ) ) == length( d.v1 ) );
                    v = v * diag( d .^ b ) / v;
                else
                    v = BaseExtDouble;
                end
            end
        end

        function v = lt( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( a2 < b2 ) );
        end

        function v = gt( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( a2 > b2 ) );
        end

        function v = le( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( a2 <= b2 ) );
        end

        function v = ge( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( a2 >= b2 ) );
        end

        function v = ne( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 ~= b1 ) | ( a2 ~= b2 );
        end

        function v = eq( a, b )
            if isa( a, 'BaseExtDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'BaseExtDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 == b1 ) & ( ~isfinite( a1 ) | ( a2 == b2 ) );
        end

        function v = colon( a, d, b )
            if nargin < 3
                b = d;
                d = 1;
            end
            if ~isa( a, 'BaseExtDouble' )
                a = a.Promote( a );
            end
            if ~isa( b, 'BaseExtDouble' )
                b = a.Promote( b );
            end
            if ~isa( d, 'BaseExtDouble' )
                d = a.Promote( d );
            end
            c = double( floor( ( b - a ) ./ d ) );
            v = a + ( 0 : c ) .* d;
        end

        function v = ctranspose( v )
            v.v1 = v.v1';
            v.v2 = v.v2';
        end

        function v = transpose( v )
            v.v1 = v.v1.';
            v.v2 = v.v2.';
        end

        function v = permute( v, DimOrder )
            v.v1 = permute( v.v1, DimOrder );
            v.v2 = permute( v.v2, DimOrder );
        end

        function v = ipermute( v, DimOrder )
            v.v1 = ipermute( v.v1, DimOrder );
            v.v2 = ipermute( v.v2, DimOrder );
        end




        function v = subsref( v, s )
            if strcmp( s(1).type, '.' )
                v = builtin( 'subsref', v, s );
            else
                v.v1 = Index( v.v1, s.subs{:} );
                v.v2 = Index( v.v2, s.subs{:} );
            end
        end

        function v = subsasgn( v, s, b )
            if strcmp( s(1).type, '.' )
                if length( s ) > 1
                    v = builtin( 'subsasgn', v, s( 1 ), subsasgn( builtin( 'subsref', v, s( 1 ) ), s( 2 : end ), b ) );
                else
                    v = builtin( 'subsasgn', v, s, b );
                end
            else
                if ~isa( v, 'BaseExtDouble' )
                    v = b.Promote( v );
                end
                if ~isa( b, 'BaseExtDouble' )
                    b = v.Promote( b );
                end
                v.v1 = Assign( v.v1, b.v1, s.subs{:} );
                v.v2 = Assign( v.v2, b.v2, s.subs{:} );
            end
        end

        function v = subsindex( v )
            v = v.v1;
        end

        function [ v, Indices ] = sort( v, varargin )
            DimIndex = find( cellfun( @isnumeric, varargin ), 1 );
            if isempty( DimIndex )
                Dim = [];
            else
                Dim = varargin{ DimIndex };
                varargin = varargin( [ 1 : ( DimIndex - 1 ), ( DimIndex + 1 ) : end ] );
            end
            CMIndex = find( strcmpi( varargin, 'ComparisonMethod' ), 1 );
            if isempty( CMIndex )
                cm = [];
            else
                cm = varargin{ CMIndex + 1 };
                varargin = varargin( [ 1 : ( CMIndex - 1 ), ( CMIndex + 2 ) : end ] );
            end
            if nargin < 2 || isempty( Dim )
                Dim = find( size( v.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            end
            if nargin < 3 || isempty( cm )
                cm = 'auto';
            end
            if isa( v, 'BaseExtDouble' )
                if strcmpi( cm( 1 ), 'r' )
                    a = real( v );
                    b = imag( v );
                elseif ( length( cm ) > 1 ) && strcmpi( cm( 1 : 2 ), 'ab' )
                    a = abs( v );
                    b = angle( v );
                else
                    if isreal( v )
                        a = real( v );
                        b = imag( v );
                    else
                        a = abs( v );
                        b = angle( v );
                    end
                end
                Size = size( v.v1 );
                if any( Size == 0 )
                    Indices = [];
                    return
                end
                Blocks = arrayfun( @( x ) ones( x, 1 ), Size, 'UniformOutput', false );
                Blocks{ Dim } = Size( Dim );
                xv1 = mat2cell( v.v1, Blocks{:} );
                xv2 = mat2cell( v.v2, Blocks{:} );
                xa1 = mat2cell( a.v1, Blocks{:} );
                xa2 = mat2cell( a.v2, Blocks{:} );
                xb1 = mat2cell( b.v1, Blocks{:} );
                xb2 = mat2cell( b.v2, Blocks{:} );
                Indices = cell( size( xv1 ) );
                for i = 1 : numel( xv1 )
                    [ ~, Indices{ i } ] = sortrows( [ xa1{ i }(:), xa2{ i }(:), xb1{ i }(:), xb2{ i }(:) ], varargin{:} );
                    xv1{ i } = Index( xv1{ i }, Indices{ i } );
                    xv2{ i } = Index( xv2{ i }, Indices{ i } );
                end
                Indices = cell2mat( Indices );
                v       = v.Make( cell2mat( xv1 ), cell2mat( xv2 ) );
            else
                if nargout > 1
                    [ v, Indices ] = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                else
                    v = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                end
                v = v.Promote( v );
            end
        end

        function v = sum( v, Dim )
            if nargin < 2
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = s.Plus( v.Make( x1{ i }, x2{ i } ) );
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = s.Plus( x{ i } );
                end
            end
            v = s;
        end

        function v = prod( v, Dim )
            if nargin < 2
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    v = ones( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = s.Times( v.Make( x1{ i }, x2{ i } ) );
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    v = ones( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = s.Times( x{ i } );
                end
            end
            v = s;
        end

        function [ s, i ] = max( a, b, Dim )
            if nargin < 3
                Dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if isempty( b )
                if isempty( a )
                    s = a.Make(0,0); s = s(zeros(0,1));
                    i = [];
                    return
                end
                if isa( a, 'BaseExtDouble' )
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a.v1 ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a.v1 );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1, Blocks{:} );
                    x2 = mat2cell( a.v2, Blocks{:} );
                    s = a.Make( x1{ 1 }, x2{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = max( a.Make( x1{ j }, x2{ j } ), s );
                        i( ii ) = j;
                    end
                else
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x = mat2cell( a, Blocks{:} );
                    s = x{ 1 };
                    for j = 2 : Length
                        s = max( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'BaseExtDouble' )
                    a = b.Promote( a );
                end
                if ~isa( b, 'BaseExtDouble' )
                    b = a.Promote( b );
                end
                [ a, b ] = BaseExtDouble.ExpandSingleton( a, b );
                i = ( a.v1 > b.v1 ) | ( ( a.v1 == b.v1 ) & ( a.v2 > b.v2 ) );
                s = b;
                s.v1 = Assign( s.v1, Index( a.v1, i ), i );
                s.v2 = Assign( s.v2, Index( a.v2, i ), i );
            end
        end

        function [ s, i ] = min( a, b, Dim )
            if nargin < 3
                Dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if isempty( b )
                if isa( a, 'BaseExtDouble' )
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a.v1 ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a.v1 );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1, Blocks{:} );
                    x2 = mat2cell( a.v2, Blocks{:} );
                    s = a.Make( x1{ 1 }, x2{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = min( a.Make( x1{ j }, x2{ j } ), s );
                        i( ii ) = j;
                    end
                else
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x = mat2cell( a, Blocks{:} );
                    s = x{ 1 };
                    for j = 2 : Length
                        s = min( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'BaseExtDouble' )
                    a = b.Promote( a );
                end
                if ~isa( b, 'BaseExtDouble' )
                    b = a.Promote( b );
                end
                [ a, b ] = BaseExtDouble.ExpandSingleton( a, b );
                i = ( a.v1 < b.v1 ) | ( ( a.v1 == b.v1 ) & ( a.v2 < b.v2 ) );
                s = b;
                s.v1 = Assign( s.v1, Index( a.v1, i ), i );
                s.v2 = Assign( s.v2, Index( a.v2, i ), i );
            end
        end

        function v = cumsum( v, Dim )
            if nargin < 2
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = s.Plus( v.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = s.Plus( x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            v = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = diff( v, Dim )
            if nargin < 2
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = v.Make( x1{ i }, x2{ i } );
                    d = t.Minus( s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = x{ i };
                    d = t.Minus( s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            end
            v = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = cumprod( v, Dim )
            if nargin < 2
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = s.Times( v.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = s.Times( x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            v = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = cummax( v, Dim )
            if nargin < 3
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = s.max( v.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = s.max( x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            v = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = cummin( v, Dim )
            if nargin < 3
                Dim = [];
            end
            if isa( v, 'BaseExtDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = s.min( v.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    v = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = s.min( x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            v = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = dot( a, b, Dim )
            if nargin < 3
                Dim = [];
            end
            if nargin < 3
                Dim = [];
            end
            if ( length( a ) == numel( a ) ) && ( length( b ) == numel( b ) )
                a = a.Vec();
                b = b.Vec();
            end
            v = a .* b;
            v = sum( v, Dim );
        end

        function v = norm( v, p )
            if nargin < 2
                p = 2;
            end
            if ( sum( size( v ) ~= 1 ) <= 1 ) || ( numel( p ) ~= 1 )
                v = abs( v );
                if ( p == 2 ) || ( numel( p ) ~= 1 )
                    v = v.Vec();
                    v = sqrt( sum( v .* v ) );
                elseif p == Inf
                    v = max( v );
                elseif p == -Inf
                    v = min( v );
                elseif p == 1
                    v = sum( v );
                else
                    v = ( sum( v .^ p ) ) .^ ( 1 ./ p );
                end
            else
                if ( p == 2 )
                    v = sqrt( max( eig( v' * v ) ) );
                elseif p == Inf
                    v = max( sum( abs( v ), 2 ) );
                elseif p == 1
                    v = max( sum( abs( v ), 1 ) );
                else
                    v = v.Promote( 0 );
                end
            end
        end

        function v = abs( v )
            if isreal( v )
                Select = v.v1 < 0;
                v.v1 = Assign( v.v1, -Index( v.v1, Select ), Select );
                v.v2 = Assign( v.v2, -Index( v.v2, Select ), Select );
            else
                real_v = real( v );
                imag_v = imag( v );
                v = sqrt( real_v .* real_v + imag_v .* imag_v );
            end
        end

        function v = sign( v )
            if isreal( v )
                v = sign( v.v1 );
            else
                abs_v = abs( v );
                [ v.v1, v.v2 ] = v.EDDividedByED( v.v1, v.v2, abs_v.v1, abs_v.v2, true );
            end
        end

        function v = floor( v )
            x1 = floor( v.v1 );
            x2 = x1 .* 0;
            Select = x1 == v.v1;
            v_sel = v.Index(Select);
            x2_sel = floor( v_sel.v2 );
            x2 = Assign( x2, x2_sel, Select );
            [ x1, x2 ] = v.Normalize( x1, x2 );
            v = v.Make( x1, x2 );
        end

        function v = ceil( v )
            x1 = ceil( v.v1 );
            x2 = x1 .* 0;
            Select = x1 == v.v1;
            v_sel = v.Index(Select);
            x2_sel = ceil( v_sel.v2 );
            x2 = Assign( x2, x2_sel, Select );
            [ x1, x2 ] = v.Normalize( x1, x2 );
            v = v.Make( x1, x2 );
        end

        function v = fix( v )
            x1 = fix( v.v1 );
            x2 = x1 .* 0;
            Select = x1 == v.v1;
            v_sel = v.Index(Select);
            x2_sel = fix( v_sel.v2 );
            x2 = Assign( x2, x2_sel, Select );
            [ x1, x2 ] = v.Normalize( x1, x2 );
            v = v.Make( x1, x2 );
        end

        function v = round( v )
            x1 = round( v.v1 );
            x2 = x1 .* 0;
            Select = x1 == v.v1;
            vSelect = v.Index(Select);
            x2Select = round( vSelect.v2 );
            x2 = Assign( x2, x2Select, Select );
            Select = ( ~Select ) & ( abs( x1 - v.v1 ) == 0.5 ) & ( v.v2 < 0 );
            x2 = Assign( x2, Index( x2, Select ) - 1, Select );
            [ x1, x2 ] = v.Normalize( x1, x2 );
            v = v.Make( x1, x2 );
        end

        function v = realsqrt( v )
            Select = v < 0;
            v.v1 = Assign( v.v1, NaN, Select );
            v.v2 = Assign( v.v2, NaN, Select );
            Select = v > 0;
            x = 1 ./ sqrt( Index( v.v1, Select ) );
            vx = Index( v.v1, Select ) .* x;
            [ t1_, t2_ ] = v.DoubleTimesUnderlying( vx, vx );
            t = v.Make( Index( v.v1, Select ), Index( v.v2, Select ) ) - v.Make( t1_, t2_ );
            [ t1_, t2_ ] = v.EDPlusUnderlying( vx, zeros( size( vx ) ), t.v1 .* ( x * 0.5 ) );
            t = v.Make( t1_, t2_ );
            v.v1 = Assign( v.v1, t.v1, Select );
            v.v2 = Assign( v.v2, t.v2, Select );
        end

        function v = sqrt( v )
            Select = v ~= 0;
            v_sel = v.Index(Select);
            x = 1 ./ sqrt( v_sel.v1 );
            vx = v_sel.v1 .* x;
            [ t1_, t2_ ] = v.DoubleTimesUnderlying( vx, vx );
            t = v_sel - v.Make( t1_, t2_ );
            [ t1_, t2_ ] = v.EDPlusUnderlying( vx, zeros( size( vx ) ), t.v1 .* ( x * 0.5 ) );
            t = v.Make( t1_, t2_ );
            v = v.Assign(t, Select);
        end

        function v = sqrtm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( sqrt( d ) ) / v;
        end

        function v = exp( v, expm1Flag )
            if nargin < 2
                expm1Flag = false;
            end
            if ~isreal( v )
                [ sin_imag_v, cos_imag_v ] = sincos( imag( v ) );
                Rotation = cos_imag_v + 1i .* sin_imag_v;
                if expm1Flag
                    v = expm1( real( v ) ) .* Rotation + Rotation;
                else
                    v = exp( real( v ) ) .* Rotation;
                end
                return
            end

            % Strategy:  We first reduce the size of x by noting that
            % exp(kr + m * log(2)) = 2^m * exp(r)^k
            % where m and k are integers.  By choosing m appropriately
            % we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
            % evaluated using the familiar Taylor series.  Reducing the
            % argument substantially speeds up the convergence.
            NSquare = v.ExpRescale;
            k = 2.0 ^ NSquare;
            inv_k = 1.0 / k;
            Threshold = inv_k .* v.tiny.v1;

            m = floor( v.v1 ./ v.log_2.v1 + 0.5 );
            r = v - v.log_2 .* m;
            r = r.TimesPowerOf2( inv_k );

            p = r .* r;
            s = r + p.TimesPowerOf2( 0.5 );
            p = p .* r;
            t = p .* v.InverseFactorial.Index( 1 );
            for i = 2 : v.NInverseFactorial
                s = s + t;
                p = p .* r;
                t = p .* v.InverseFactorial.Index( i );
                if all( abs( Index( t.v1, ':' ) ) <= Threshold )
                    break
                end
            end

            s = s + t;

            for i = 1 : NSquare
                s = s.TimesPowerOf2( 2.0 ) + s .* s;
            end

            if expm1Flag
                Select = m ~= 0;
                [ t1, t2 ] = v.EDPlusUnderlying( Index( s.v1, Select ), Index( s.v2, Select ), 1.0 );
                s.v1 = Assign( s.v1, t1, Select );
                s.v2 = Assign( s.v2, t2, Select );
            else
                s = s + 1.0;
            end

            v = v.Make( pow2( s.v1, m ), pow2( s.v2, m ) );
            if expm1Flag
                [ t1, t2 ] = v.EDPlusUnderlying( Index( v.v1, Select ), Index( v.v2, Select ), -1.0 );
                v.v1 = Assign( v.v1, t1, Select );
                v.v2 = Assign( v.v2, t2, Select );
            end
        end

        function v = expm1( v )
            v = exp( v, true );
        end

        function v = expm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( exp( d ) ) / v;
        end

        function x = log( v )
            x = v.Make( log( v.v1 ), zeros( size( v.v1 ) ) );
            for i = 1 : v.LogSteps
                x = x + v .* exp( -x ) - 1.0;
            end
        end

        function v = log2( v )
            v = log( v ) ./ v.log_2;
        end

        function v = log10( v )
            v = log( v ) ./ v.log_10;
        end

        function v = logm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( log( d ) ) / v;
        end

        function v = funm( v, f )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( f( d ) ) / v;
        end

        function [ sin_v, cos_v ] = sincos( v )
            if ~isreal( v )
                exp_Piv = exp( v.TimesPowerOf2( 1i ) );
                exp_Niv = 1 ./ exp_Piv;
                sin_v = exp_Piv - exp_Niv;
                sin_v = sin_v.TimesPowerOf2( -0.5i );
                cos_v = exp_Piv + exp_Niv;
                cos_v = cos_v.TimesPowerOf2( +0.5 );
                return
            end
            % Strategy.  To compute sin(x), cos(x), we choose integers a, b so that
            % x = s + a * (pi/2) + b * (pi/16)
            % and |s| <= pi/32.  Using the fact that
            % sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
            % we can compute sin(x) from sin(s), cos(s).  This greatly increases the convergence of the sine Taylor series.

            z = round( v ./ v.piT2 );
            r = v - v.piT2 .* z;

            q = floor( r.v1 ./ v.piD2.v1 + 0.5 );
            t = r - v.piD2 .* q;
            j = q;
            abs_j = abs( j );

            q = floor( t.v1 ./ v.piD16.v1 + 0.5 );
            t = t - v.piD16 .* q;
            k = q;
            abs_k = abs( k );

            test = ( j >= -2 ) & ( j <= 2 );
            assert( all( test(:) ) );
            test = abs_k <= 4;
            assert( all( test(:) ) );

            [ sin_t, cos_t ] = t.SinCosTaylor();

            sin_v = sin_t;
            cos_v = cos_t;

            a = v.CosTable.Index( double(abs_k + 1) );
            b = v.SinTable.Index( double(abs_k + 1) );

            a = reshape( a, size( v ) );
            b = reshape( b, size( v ) );

            a_sin_t = a .* sin_t;
            b_sin_t = b .* sin_t;
            a_cos_t = a .* cos_t;
            b_cos_t = b .* cos_t;

            Select = k > 0;
            if any(Select(:))
                a_sin_t_s = a_sin_t.Index(Select);
                a_cos_t_s = a_cos_t.Index(Select);
                b_sin_t_s = b_sin_t.Index(Select);
                b_cos_t_s = b_cos_t.Index(Select);
                sin_v = sin_v.Assign( a_sin_t_s + b_cos_t_s, Select );
                cos_v = cos_v.Assign( a_cos_t_s - b_sin_t_s, Select );
            end

            Select = k < 0;
            if any(Select(:))
                a_sin_t_s = a_sin_t.Index(Select);
                a_cos_t_s = a_cos_t.Index(Select);
                b_sin_t_s = b_sin_t.Index(Select);
                b_cos_t_s = b_cos_t.Index(Select);
                sin_v = sin_v.Assign( a_sin_t_s - b_cos_t_s, Select );
                cos_v = cos_v.Assign( a_cos_t_s + b_sin_t_s, Select );
            end

            Select = j == 1;
            if any(Select(:))
                cos_v_s = cos_v.Index(Select);
                sin_v_s = sin_v.Index(Select);
                sin_v = sin_v.Assign( cos_v_s, Select );
                cos_v = cos_v.Assign( -sin_v_s, Select );
            end

            Select = j == -1;
            if any(Select(:))
                cos_v_s = cos_v.Index(Select);
                sin_v_s = sin_v.Index(Select);
                sin_v = sin_v.Assign( -cos_v_s, Select );
                cos_v = cos_v.Assign( sin_v_s, Select );
            end

            Select = abs_j == 2;
            if any(Select(:))
                sin_v_s = sin_v.Index(Select);
                cos_v_s = cos_v.Index(Select);
                sin_v = sin_v.Assign( -sin_v_s, Select );
                cos_v = cos_v.Assign( -cos_v_s, Select );
            end
        end

        function v = sin( v )
            [ v, ~ ] = sincos( v );
        end

        function v = asin( v )
            assert( all( abs( Index( v.v1, ':' ) ) <= 1 ) );
            v = atan2( v, sqrt( max( 1 - v.*v, 0 ) ) );
        end

        function v = cos( v )
            [ ~, v ] = sincos( v );
        end

        function v = acos( v )
            assert( all( abs( Index( v.v1, ':' ) ) <= 1 ) );
            v = atan2( sqrt( max( 1 - v.*v, 0 ) ), v );
        end

        function v = tan( v )
            [ sin_v, cos_v ] = sincos( v );
            v = sin_v ./ cos_v;
        end

        function v = atan( v )
            v = atan2( v, v.Make( 1, 0 ) );
        end

        function v = atan2( y, x )
            r = sqrt( x.*x + y.*y );
            xx = x ./ r;
            yy = y ./ r;
            Select = abs( xx.v1 ) > abs( yy.v1 );
            v = y.Promote( atan2( y.v1, x.v1 ) );
            [ sin_z, cos_z ] = sincos( v );
            t = yy;

            if any(Select(:))
                tSelect = t.Index(Select);
                sinZSelect = sin_z.Index(Select);
                cosZSelect = cos_z.Index(Select);
                t = t.Assign( ( tSelect - sinZSelect ) ./ cosZSelect, Select );
            end

            Select = ~Select;
            if any(Select(:))
                xxSelect = xx.Index(Select);
                sinZSelect = sin_z.Index(Select);
                cosZSelect = cos_z.Index(Select);
                t = t.Assign( ( xxSelect - cosZSelect ) ./ ( -sinZSelect ), Select );
            end
            v = v + t;
        end

        function v = sinh( v )
            exp_v = exp( v );
            v = exp_v - 1 ./ exp_v;
            v = v.TimesPowerOf2( 0.5 );
        end

        function v = asinh( v )
            v = log( v + sqrt( v.*v + 1 ) );
        end

        function v = cosh( v )
            exp_v = exp( v );
            v = exp_v + 1 ./ exp_v;
            v = v.TimesPowerOf2( 0.5 );
        end

        function v = acosh( v )
            v = log( v + sqrt( v.*v - 1 ) );
        end

        function [ sinh_v, cosh_v ] = sinhcosh( v )
            exp_Pv = exp( v );
            exp_Nv = 1 ./ exp_Pv;
            sinh_v = exp_Pv - exp_Nv;
            sinh_v = sinh_v.TimesPowerOf2( 0.5 );
            cosh_v = exp_Pv + exp_Nv;
            cosh_v = cosh_v.TimesPowerOf2( 0.5 );
        end

        function v = tanh( v )
            [ sinh_v, cosh_v ] = sinhcosh( v );
            v = sinh_v ./ cosh_v;
        end

        function v = atanh( v )
            v = log( ( 1 + v ) ./ ( 1 - v ) );
            v = v.TimesPowerOf2( 0.5 );
        end

        function v = mod( v, b )
            v = v - b .* floor( v ./ b );
        end

        function v = rem( v, b )
            v = v - b .* fix( v ./ b );
        end

        function [ v, U, p ] = lu( v, type )
            if nargin < 2
                [ v, U, p ] = luExt( v );
            else
                [ v, U, p ] = luExt( v, type );
            end
        end

        function [ q, v ] = qr( v )
            [ q, v ] = qrExt( v );
        end

        function v = det( v )
            [ m, n ] = size( v );
            if m ~= n
                throw( MException( 'MATLAB:square', 'Matrix must be square.' ) );
            end
            [ ~, u, P ] = lu( v );
            DetP = det( P );
            if DetP > 0
                v = prod( diag( u ) );
            elseif DetP < 0
                v = -prod( diag( u ) );
            else
                v = v.Make( NaN, NaN ); % scalar NaN, can't use 'like' here
            end
        end

        function v = inv( v )
            n = size( v, 1 );
            v = v \ v.Make( eye( n ), zeros( n ) );
        end

        function [ v, p ] = chol( v, type )
            [ v, d ] = ldl( v, 'vector_d' );
            v = v .* sqrt( d.' );
            if any( d < 0 )
                p = 1;
            else
                p = 0;
            end
            if nargin < 2 || strcmp( type, 'upper' )
                v = v.';
            end
        end

        function [ L, d ] = ldl( v, type )
            if nargin < 2
                type = [];
            end
            [ L, d ] = ldlExt( v, type );
        end

        function [ v, d ] = eig( x )
            [ v, d ] = eigExt( x );
        end

        function w = conv( u, v )
            w = convExt( u, v );
        end

        function [ C, ia, ic ] = unique( A, varargin )
            Rows = strcmpi( varargin, 'rows');
            if any( Rows )
                varargin( Rows ) = [];
                RowFlag = 0;
            else
                RowFlag = size( A, 1 ) == 1;
                A = A.Vec();
            end
            Size = size( A );
            A_orig = A;
            A = [ A.v1, A.v2 ];
            [ C, ia, ic ] = unique( A, 'rows', varargin{:} );
            n = Size( 2 );
            C = A_orig.Make( Index( C, ':', 1:n ), Index( C, ':', ( n + 1 ) : ( 2 * n ) ) );
            if RowFlag
                C = C.';
            end
        end

        function v = mean( v, Dim )
            if ( nargin < 2 ) || isempty( Dim )
                Dim = find( size( v.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            elseif strcmpi( Dim, 'all' )
                v = v.Vec();
                Dim = 1;
            end
            Size = size( v.v1 );
            n = prod( Size( Dim ) );
            v = sum( v, Dim ) ./ n;
        end

        function v = median( v, Dim )
            if ( nargin < 2 ) || isempty( Dim )
                Dim = find( size( v.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            elseif strcmpi( Dim, 'all' )
                v = v.Vec();
                Dim = 1;
            end
            Size = size( v.v1 );
            n = prod( Size( Dim ) );

            if n == 0
                Size( Dim ) = 0;
                v = NaN( Size, 'like', v );
                return;
            end

            NotDim = setdiff( 1 : numel( Size ), Dim );
            v = reshape( permute( v, [ Dim, NotDim ] ), [ n, Size( NotDim ) ] );

            v = sort( v, 1 );

            if mod( n, 2 ) == 1
                Middle = ( n + 1 ) * 0.5;
                v = v.Make( Index( v.v1, Middle, ':' ), Index( v.v2, Middle, ':' ) );
            else
                Middle1 = n * 0.5;
                Middle2 = Middle1 + 1;
                m1 = v.Make( Index( v.v1, Middle1, ':' ), Index( v.v2, Middle1, ':' ) );
                m2 = v.Make( Index( v.v1, Middle2, ':' ), Index( v.v2, Middle2, ':' ) );
                v = 0.5 * ( m1 + m2 );
            end
            v = ipermute( reshape( v, [ ones( 1, numel( Dim ) ), Size( NotDim ) ] ), [ Dim, NotDim ] );
        end

        function v = std( v, varargin )
            v = sqrt( var( v, varargin{:} ) );
        end

        function v = var( v, Flag, Dim )
            if nargin < 3
                Dim = [];
                if nargin < 2
                    Flag = 0;
                end
            end

            if isempty( Dim )
                Dim = find( size( v.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            end

            Size = size( v.v1 );
            n = prod( Size( Dim ) );

            if ( n == 0 ) || ( ( n == 1 ) && ( Flag == 0 ) )
                Size( Dim ) = 1;
                v = NaN( Size, 'like', v );
                return;
            end

            Mu = mean( v, Dim );
            v = v - Mu;
            v = v .* v;
            v = sum( v, Dim );

            if Flag == 1
                v = v ./ n;
            else
                v = v ./ ( n - 1 );
            end
        end

        function [ X, Y ] = meshgrid( x, y )
            [ X1, Y1 ] = meshgrid( x.v1, y.v1 );
            [ X2, Y2 ] = meshgrid( x.v2, y.v2 );
            X = x.Make( X1, X2 );
            Y = x.Make( Y1, Y2 );
        end

        function y = linspace( a, b, n )
            if nargin < 3
                n = 100;
            end

            if ~isa( a, 'BaseExtDouble' )
                a = b.Promote( a );
            end

            if ~isa( b, 'BaseExtDouble' )
                b = a.Promote( b );
            end

            if n < 1
                y = zeros( 0, 1, 'like', a );
                return;
            end

            if n == 1
                y = b;
                return;
            end

            Step = ( b - a ) ./ ( n - 1 );
            Indices = 0 : ( n - 1 );
            y = a + Step .* Indices;

            if n > 1
                y.v1 = Assign( y.v1, b.v1, length(y.v1) );
                y.v2 = Assign( y.v2, b.v2, length(y.v2) );
            end
        end

        function v = MTimes( a, b )
            [ R, c ] = size( a );
            [ r, C ] = size( b );
            if ( ( R == 1 ) && ( c == 1 ) ) || ( ( r == 1 ) && ( C == 1 ) )
                v = a.Times( b );
                return
            end
            v = zeros( R, C, 'like', a );
            if isa( b, 'BaseExtDouble' )
                for c = 1 : C
                    t = sum( a .* b.Make( Index( b.v1, ':', c ).', Index( b.v2, ':', c ).' ), 2 );
                    v.v1 = Assign( v.v1, t.v1, ':', c );
                    v.v2 = Assign( v.v2, t.v2, ':', c );
                end
            else
                for c = 1 : C
                    t = sum( a .* b( :, c ).', 2 );
                    v.v1 = Assign( v.v1, t.v1, ':', c );
                    v.v2 = Assign( v.v2, t.v2, ':', c );
                end
            end
        end

        function v = MLDivide( a, v )
            [ Ra, Ca ] = size( a );
            [ Rv, Cv ] = size( v );
            if ( ( Ra == 1 ) && ( Ca == 1 ) ) || ( ( Rv == 1 ) && ( Cv == 1 ) )
                v = a.LDivide( v );
                return
            end
            if ~isa( v, 'BaseExtDouble' )
                v = BaseExtDouble( v );
            end
            assert( Ra == Rv );
            if Ra ~= Ca
                [ q, a ] = qr( a );
                v = q' * v;
                v = BackSubstitutionExt( v, a );
                return
            end
            if BaseExtDouble.IsEqualWithExpansion( triu( a, 1 ), 0 )
                % Lower triangular
                v = ForwardEliminationExt( v, a );
                return
            elseif BaseExtDouble.IsEqualWithExpansion( tril( a, -1 ), 0 )
                % Upper triangular
                v = BackSubstitutionExt( v, a );
                return
            elseif BaseExtDouble.IsEqualWithExpansion( a, a' )
                [ L, d ] = ldl( a, 'vector_d' );
                if all( all( isfinite( L ) ) ) && all( isfinite( d ) )
                    % Positive definite
                    v = ForwardEliminationExt( v, L );
                    v = v ./ d;
                    v = BackSubstitutionExt( v, L' );
                    return
                end
            end
            % Triangular factorization
            [ L, U, p ] = lu( a, 'vector' );

            % Permutation and forward elimination
            v.v1 = Index( v.v1, p, ':' );
            v.v2 = Index( v.v2, p, ':' );
            v = ForwardEliminationExt( v, L );

            % Back substitution
            v = BackSubstitutionExt( v, U );
        end

        function v = MRDivide( v, a )
            v = MLDivide(a', v' )';
        end

        function [ v, Indices ] = Sort( v, Dim, cm, varargin )
            if nargin < 2 || isempty( Dim )
                Dim = find( size( v.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            end
            if nargin < 3 || isempty( cm )
                cm = 'auto';
            end
            if isa( v, 'BaseExtDouble' )
                if strcmpi( cm( 1 ), 'r' )
                    a = real( v );
                    b = imag( v );
                elseif ( length( cm ) > 1 ) && strcmpi( cm( 1 : 2 ), 'ab' )
                    a = abs( v );
                    b = angle( v );
                else
                    if isreal( v )
                        a = real( v );
                        b = imag( v );
                    else
                        a = abs( v );
                        b = angle( v );
                    end
                end
                Size = size( v.v1 );
                if any( Size == 0 )
                    Indices = [];
                    return
                end
                Blocks = arrayfun( @( x ) ones( x, 1 ), Size, 'UniformOutput', false );
                Blocks{ Dim } = Size( Dim );
                xv1 = mat2cell( v.v1, Blocks{:} );
                xv2 = mat2cell( v.v2, Blocks{:} );
                xa1 = mat2cell( a.v1, Blocks{:} );
                xa2 = mat2cell( a.v2, Blocks{:} );
                xb1 = mat2cell( b.v1, Blocks{:} );
                xb2 = mat2cell( b.v2, Blocks{:} );
                Indices = cell( size( xv1 ) );
                for i = 1 : numel( xv1 )
                    [ ~, Indices{ i } ] = sortrows( [ xa1{ i }(:), xa2{ i }(:), xb1{ i }(:), xb2{ i }(:) ], varargin{:} );
                    xv1{ i } = xv1{ i }( Indices{ i } );
                    xv2{ i } = xv2{ i }( Indices{ i } );
                end
                Indices = cell2mat( Indices );
                v       = v.Make( cell2mat( xv1 ), cell2mat( xv2 ) );
            else
                if nargout > 1
                    [ v, Indices ] = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                else
                    v = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                end
                v = BaseExtDouble( v );
            end
        end



        function c = Diff( v, Dim )
            if isa( v, 'BaseExtDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = v.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = v.Make( x1{ i }, x2{ i } );
                    d = t.Minus( s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            else
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( Dim );
                if Length == 0
                    c = zeros( Size, 'like', v );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = x{ i };
                    d = t.Minus( s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            end
            c = v.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = Norm( v, p )
            if ( sum( size( v ) ~= 1 ) <= 1 ) || ( numel( p ) ~= 1 )
                v = abs( v );
                if ( nargin < 2 ) || ( p == 2 ) || ( numel( p ) ~= 1 )
                    v = Vec( v );
                    v = sqrt( sum( v .* v ) );
                elseif p == Inf
                    v = max( v );
                elseif p == -Inf
                    v = min( v );
                elseif p == 1
                    v = sum( v );
                else
                    v = ( sum( v .^ p ) ) .^ ( 1 ./ p );
                end
            else
                if ( nargin < 2 ) || ( p == 2 )
                    v = sqrt( max( eig( v' * v ) ) );
                elseif p == Inf
                    v = max( sum( abs( v ), 2 ) );
                elseif p == 1
                    v = max( sum( abs( v ), 1 ) );
                else
                    v = BaseExtDouble;
                end
            end
        end

        function v = Vec( v )
            v.v1 = Vec( v.v1 );
            v.v2 = Vec( v.v2 );
        end


        function v = pow2( v, b )
            if nargin == 1
                v = pow2( 2, v );
            else
                if isa( v, 'BaseExtDouble' )
                    v.v1 = pow2( v.v1, double( b ) );
                    v.v2 = pow2( v.v2, double( b ) );
                else
                    v = pow2( v, double( b ) );
                end
            end
        end

        function v = TimesPowerOf2( v, b )
            if isa( v, 'BaseExtDouble' )
                v.v1 = TimesPowerOf2( v.v1, b );
                v.v2 = TimesPowerOf2( v.v2, b );
            else
                v = TimesPowerOf2( v, double( b ) );
            end
        end

        function v = SinTaylor( v )
            Threshold = 0.5 .* abs( Index( v.v1, ':' ) ) .* v.tiny.v1;
            x = - v .* v;
            r = v;
            for i = 1 : 2 : v.NInverseFactorial
                r = r .* x;
                t = r .* v.InverseFactorial.Index( i );
                v = v + t;
                if all( abs( Index( t.v1, ':' ) ) <= Threshold )
                    break
                end
            end
        end

        function [ sin_v, cos_v ] = SinCosTaylor( v )
            sin_v = SinTaylor( v );
            cos_v = sqrt( 1 - sin_v .* sin_v );
        end

        function [ a1, a2 ] = Split( a )
            [ c1, c2 ] = Split( a.v1 );
            [ c3, c4 ] = Split( a.v2 );
            a1 = a.Make( c1, c2 );
            a2 = a.Make( c3, c4 );
        end

        function v = ToRand( v ) % Fills v with random values in place.
            Size = size( v.v1 );
            v.v1 = rand( Size, 'like', v.v1 );
            v.v2 = v.tiny.v1 .* ( rand( Size, 'like', v.v1 ) - 0.5 );
        end

        function v = ToRandn( v ) % Fills v with normal random values in place.
            Size = size( v.v1 );
            v.v1 = Vec( v.v1 );
            v.v2 = Vec( v.v2 );
            if mod( numel( v ), 2 ) ~= 0
                v = [ v; v.Promote( 0 ) ];
                Expanded = true;
            else
                Expanded = false;
            end
            N = numel( v );
            H = 0.5 * N;
            v = reshape( v, H, 2 );
            v = ToRand( v );
            U = Index( v, ':', 1 );
            V = Index( v, ':', 2 );
            R = sqrt( -2 * log( U ) );
            Theta = 2 * v.pi * V;
            [ S, C ] = sincos( Theta );
            v = Assign( v, R .* C, ':', 1 );
            v = Assign( v, R .* S, ':', 2 );
            if Expanded
                v = Assign( Vec( v ), [], N );
            end
            v = reshape( v, Size );
        end

        function v = LDivide( b, a )
            v = RDivide( a, b );
        end

        function v = Minus( a, b )
            v = Plus( a, -b );
        end

    end

    methods % Not sealed

        function v = Plus( a, b )
            if isa( a, 'BaseExtDouble' )
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDPlusED( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = BaseExtDouble.EDPlusUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                end
                v = a.Make( x1, x2 );
            else
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDPlusUnderlying( b.v1, b.v2, Promote( b.v1, a ) );
                else
                    % [ x1, x2 ] = BaseExtDouble.DoublePlusUnderlying( double( a ), double( b ) );
                    error( 'Should never happen.' );
                end
                v = b.Make( x1, x2 );
            end
        end

        function v = Times( a, b )
            if isa( a, 'BaseExtDouble' )
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDTimesED( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = BaseExtDouble.EDTimesUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                end
                v = a.Make( x1, x2 );
            else
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDTimesUnderlying( b.v1, b.v2, Promote( b.v1, a ) );
                else
                    % [ x1, x2 ] = BaseExtDouble.DoubleTimesUnderlying( double( a ), double( b ) );
                    error( 'Should never happen.' );
                end
                v = b.Make( x1, x2 );
            end
        end

        function v = RDivide( a, b )
            if isa( a, 'BaseExtDouble' )
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDDividedByED( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = BaseExtDouble.EDDividedByUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                end
                v = a.Make( x1, x2 );
            else
                if isa( b, 'BaseExtDouble' )
                    [ x1, x2 ] = BaseExtDouble.EDDividedByED( Promote( b.v1, a ), Promote( b.v1, zeros( size( a ) ) ), b.v1, b.v2 );
                else
                    % [ x1, x2 ] = BaseExtDouble.DoubleDividedByUnderlying( double( a ), double( b ) );
                    error( 'Should never happen.' );
                end
                v = b.Make( x1, x2 );
            end
        end


    end

    methods ( Static, Access = private )


        function [ s1, s2 ] = Normalize( a1, a2 )
            s1 = a1 + a2;
            t = s1 - a1;
            s2 = a2 - t;
        end

        function [ s1, s2 ] = EDPlusED( a1, a2, b1, b2 )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = BaseExtDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            [ s1, s2 ] = BaseExtDouble.DoublePlusUnderlying( a1, b1 );
            [ t1, t2 ] = BaseExtDouble.DoublePlusUnderlying( a2, b2 );
            s2 = s2 + t1;
            [ s1, s2 ] = BaseExtDouble.Normalize( s1, s2 );
            s2 = s2 + t2;
            [ s1, s2 ] = BaseExtDouble.Normalize( s1, s2 );
        end

        function [ s1, s2 ] = EDPlusUnderlying( a1, a2, b )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = BaseExtDouble.ExpandSingleton( a1, a2, b );
            end
            [ s1, s2 ] = BaseExtDouble.DoublePlusUnderlying( a1, b );
            s2 = s2 + a2;
            [ s1, s2 ] = BaseExtDouble.Normalize( s1, s2 );
        end

        function [ s1, s2 ] = DoublePlusUnderlying( a, b )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a, b ] = BaseExtDouble.ExpandSingleton( a, b );
            end
            s1 = a + b;
            bb = s1 - a;
            t11 = s1 - bb;
            t2 = b - bb;
            t1 = a - t11;
            s2 = t1 + t2;
        end

        function [ p1, p2 ] = EDTimesED( a1, a2, b1, b2 )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = BaseExtDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            [ p1, p2 ] = BaseExtDouble.DoubleTimesUnderlying( a1, b1 );
            t = a1 .* b2 + a2 .* b1;
            p2 = p2 + t;
            [ p1, p2 ] = BaseExtDouble.Normalize( p1, p2 );
        end

        function [ p1, p2 ] = EDTimesUnderlying( a1, a2, b )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = BaseExtDouble.ExpandSingleton( a1, a2, b );
            end
            [ p1, p2 ] = BaseExtDouble.DoubleTimesUnderlying( a1, b );
            p2 = p2 + a2 .* b;
            [ p1, p2 ] = BaseExtDouble.Normalize( p1, p2 );
        end

        function [ p1, p2 ] = DoubleTimesUnderlying( a, b )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a, b ] = BaseExtDouble.ExpandSingleton( a, b );
            end
            p1 = a .* b;
            [ a1, a2 ] = Split( a );
            [ b1, b2 ] = Split( b );
            t1 = a1 .* b1 - p1;
            t2 = t1 + a1 .* b2 + a2 .* b1;
            p2 = t2 + a2 .* b2;
        end

        function [ r1, r2 ] = EDDividedByED( a1, a2, b1, b2, AnySolutionWillDo )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = BaseExtDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            if nargin < 5
                AnySolutionWillDo = false;
            end
            q1 = a1 ./ b1;
            [ p1, p2 ] = BaseExtDouble.EDTimesUnderlying( b1, b2, q1 );
            [ r1, r2 ] = BaseExtDouble.EDPlusED( a1, a2, -p1, -p2 );
            q2 = r1 ./ b1;
            [ p1, p2 ] = BaseExtDouble.EDTimesUnderlying( b1, b2, q2 );
            [ r1, ~  ] = BaseExtDouble.EDPlusED( r1, r2, -p1, -p2 );
            q3 = r1 ./ b1;
            [ q1, q2 ] = BaseExtDouble.Normalize( q1, q2 );
            [ r1, r2 ] = BaseExtDouble.EDPlusED( q1, q2, q3, zeros( size( q3 ) ) );
            Select = ( b1 == 0 ) & ( b2 == 0 );
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
                    Select = repmat( Select, size( a1 ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
                    a1 = repmat( a1, size( Select ) );
                end
                a1Select = a1( Select );
                a1SelectSelect = a1Select == 0;
                a1Select = sign( a1Select ) .* Inf;
                if AnySolutionWillDo
                    a1Select = Assign( a1Select, 0, a1SelectSelect );
                end
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
            Select = isinf( b1 );
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
                    Select = repmat( Select, size( a1 ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
                    a1 = repmat( a1, size( Select ) );
                end
                a1Select = a1( Select );
                a1SelectSelect = ~isfinite( a1Select );
                a1Select = 0;
                a1Select = Assign( a1Select, NaN, a1SelectSelect );
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
        end

        function [ r1, r2 ] = EDDividedByUnderlying( a1, a2, b, AnySolutionWillDo )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = BaseExtDouble.ExpandSingleton( a1, a2, b );
            end
            if nargin < 4
                AnySolutionWillDo = false;
            end
            r1 = a1 ./ b;
            [ p1, p2 ] = BaseExtDouble.DoubleTimesUnderlying( r1, b );
            [ s, e ] = BaseExtDouble.DoublePlusUnderlying( a1, -p1 );
            e = e + a2;
            e = e - p2;
            t = s + e;
            r2 = t ./ b;
            [ r1, r2 ] = BaseExtDouble.Normalize( r1, r2 );
            Select = b == 0;
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
                    Select = repmat( Select, size( a1 ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
                    a1 = repmat( a1, size( Select ) );
                end
                a1Select = a1( Select );
                a1SelectSelect = a1Select == 0;
                a1Select = sign( a1Select ) .* Inf;
                if AnySolutionWillDo
                    a1Select = Assign( a1Select, 0, a1SelectSelect );
                end
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
            Select = isinf( b );
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a1 ) > 1 )
                    Select = repmat( Select, size( a1 ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a1 ) )
                    a1 = repmat( a1, size( Select ) );
                end
                a1Select = a1( Select );
                a1SelectSelect = ~isfinite( a1Select );
                a1Select = 0;
                a1Select = Assign( a1Select, NaN, a1SelectSelect );
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
        end

        function [ r1, r2 ] = DoubleDividedByUnderlying( a, b, AnySolutionWillDo )
            if BaseExtDouble.SingletonExpansionNotSupported
                [ a, b ] = BaseExtDouble.ExpandSingleton( a, b );
            end
            if nargin < 3
                AnySolutionWillDo = false;
            end
            r1 = a ./ b;
            [ p1, p2 ] = BaseExtDouble.DoubleTimesUnderlying( r1, b );
            [ s, e ] = BaseExtDouble.DoublePlusUnderlying( a, -p1 );
            e = e - p2;
            t = s + e;
            r2 = t ./ b;
            [ r1, r2 ] = BaseExtDouble.Normalize( r1, r2 );
            Select = b == 0;
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a ) > 1 )
                    Select = repmat( Select, size( a ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a ) )
                    a = repmat( a, size( Select ) );
                end
                a1Select = a( Select );
                a1SelectSelect = a1Select == 0;
                a1Select = sign( a1Select ) .* Inf;
                if AnySolutionWillDo
                    a1Select = Assign( a1Select, 0, a1SelectSelect );
                end
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
            Select = isinf( b );
            if any( Select(:) )
                if ( isscalar( Select ) ) && ( numel( a ) > 1 )
                    Select = repmat( Select, size( a ) );
                elseif ( numel( Select ) > 1 ) && ( isscalar( a ) )
                    a = repmat( a, size( Select ) );
                end
                a1Select = a( Select );
                a1SelectSelect = ~isfinite( a1Select );
                a1Select = 0;
                a1Select = Assign( a1Select, NaN, a1SelectSelect );
                r1 = Assign( r1, a1Select, Select );
                r2 = Assign( r2, a1Select, Select );
            end
        end

    end


    methods ( Static, Access = public )

        function Supported = TestSingletonExpansion()
            try
                a = [ 1, 1 ];
                b = [ 1; 1 ];
                c = a + b;
                Supported = all( [ size( b ), size( c ) ] == 2 );
            catch
                Supported = false;
            end
        end

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

        function v = IsEqualWithExpansion( a, b, varargin )
            v = a == b;
            v = all( v(:) );
            if nargin > 2
                for i = 1 : length( varargin )
                    if ~v
                        break
                    end
                    v = v && BaseExtDouble.IsEqualWithExpansion( a, varargin{i} );
                end
            end
        end
    end

end
