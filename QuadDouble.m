% The below code is derived from the QD C++ library.

% QD is Copyright ( c ) 2003-2009, The Regents of the University of California, through Lawrence Berkeley National Laboratory ( subject to receipt of any required approvals from U.S. Dept. of Energy ) All rights reserved.

% QD is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% ( 3 ) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 3. You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ( "Enhancements" ) to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.

% The implementation of the LUP decomposition and the backslack operator here is derived from Cleve Moler's code from "Numerical Computing with MATLAB".

% "Numerical Computing with MATLAB" is Copyright ( c ) 2004, Cleve Moler and Copyright ( c ) 2016, The MathWorks, Inc.

% "Numerical Computing with MATLAB" is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% ( 3 ) In all cases, the software is, and all modifications and derivatives of the software shall be, licensed to you solely for use in conjunction with MathWorks products and service offerings.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% The implementation of the LDL decomposition here is derived from Brian Borchers's ldlt library.

% Brian Borchers's ldlt library is Copyright ( c ) 2009, Brian Borchers.

% Brian Borchers's ldlt library is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% ( 1 ) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% ( 2 ) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
% ( 3 ) Neither the name of the New Mexico Inst of Mining & Tech nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION ) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ( INCLUDING NEGLIGENCE OR OTHERWISE ) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% WARNING: This code was primarily written by AI, based on my DoubleDouble class. Having not gone through it line by line, I cannot 100% guarantee its correctness.

classdef QuadDouble < QuadDoubleSlow

    methods
        function v = QuadDouble( in, varargin )
            if nargin == 0
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = QuadDouble.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDouble' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            elseif isa( in, 'BaseExtDouble' )
                if isprop(in, 'v3')
                    C = cell( 1, 4 );
                else
                    C = cell( 1, 2 );
                end
                [ C{ : } ] = ToSumOfDoubles( in );
                if length( C ) < 2; C{ 2 } = zeros( size( C{ 1 } ) ); end
                if length( C ) < 3; C{ 3 } = zeros( size( C{ 1 } ) ); end
                if length( C ) < 4; C{ 4 } = zeros( size( C{ 1 } ) ); end
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );
            else
                v.v1 = DoubleDouble.MakeStatic( double( in ), zeros( size( in ) ) );
                v.v2 = DoubleDouble.MakeStatic( zeros( size( in ) ), zeros( size( in ) ) );
            end
        end

        function disp( v )
            if isempty( v.v1.v1 )
                disp( '     []' );
            else
                disp( v.v1.v1 );
            end
            disp( '     +' );
            if isempty( v.v1.v2 )
                disp( '     []' );
            else
                disp( v.v1.v2 );
            end
            disp( '     +' );
            if isempty( v.v2.v1 )
                disp( '     []' );
            else
                disp( v.v2.v1 );
            end
            disp( '     +' );
            if isempty( v.v2.v2 )
                disp( '     []' );
            else
                disp( v.v2.v2 );
            end
            disp( ' ' );
        end

        function [ x1, x2, x3, x4 ] = ToSumOfDoubles( v )
            x1 = v.v1.v1;
            x2 = v.v1.v2;
            x3 = v.v2.v1;
            x4 = v.v2.v2;
        end

        function v = double( v )
            v = v.v1.v1;
        end

        function v = isscalar( v )
            v = isscalar( v.v1.v1 );
        end

        function v = isreal( v )
            v = isreal( v.v1.v1 ) && isreal( v.v1.v2 ) && isreal( v.v2.v1 ) && isreal( v.v2.v2 );
        end

        function v = isnumeric( v )
            v = isnumeric( v.v1.v1 ) && isnumeric( v.v1.v2 ) && isnumeric( v.v2.v1 ) && isnumeric( v.v2.v2 );
        end

        function v = isfinite( v )
            v = isfinite( v.v1.v1 ) & isfinite( v.v1.v2 ) & isfinite( v.v2.v1 ) & isfinite( v.v2.v2 );
        end

        function v = isinf( v )
            v = isinf( v.v1.v1 ) | isinf( v.v1.v2 ) | isinf( v.v2.v1 ) | isinf( v.v2.v2 );
        end

        function v = isnan( v )
            v = isnan( v.v1.v1 ) | isnan( v.v1.v2 ) | isnan( v.v2.v1 ) | isnan( v.v2.v2 );
        end

        function v = sparse( i, j, v, m, n, nz )
            if nargin == 0
                v = QuadDouble.Make( sparse( [] ), sparse( [] ), sparse( [] ), sparse( [] ) );
            elseif nargin < 3
                assert( nargin == 1 );
                v = i;
                v.v1.v1 = sparse( v.v1.v1 );
                v.v1.v2 = sparse( v.v1.v2 );
                v.v2.v1 = sparse( v.v2.v1 );
                v.v2.v2 = sparse( v.v2.v2 );
            elseif nargin == 3
                v.v1.v1 = sparse( i, j, v.v1.v1 );
                v.v1.v2 = sparse( i, j, v.v1.v2 );
                v.v2.v1 = sparse( i, j, v.v2.v1 );
                v.v2.v2 = sparse( i, j, v.v2.v2 );
            elseif nargin == 4
                v.v1.v1 = sparse( i, j, v.v1.v1, m, n );
                v.v1.v2 = sparse( i, j, v.v1.v2, m, n );
                v.v2.v1 = sparse( i, j, v.v2.v1, m, n );
                v.v2.v2 = sparse( i, j, v.v2.v2, m, n );
            else
                v.v1.v1 = sparse( i, j, v.v1.v1, m, n, nz );
                v.v1.v2 = sparse( i, j, v.v1.v2, m, n, nz );
                v.v2.v1 = sparse( i, j, v.v2.v1, m, n, nz );
                v.v2.v2 = sparse( i, j, v.v2.v2, m, n, nz );
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
                    varargout{ 3 } = QuadDouble.Make( v.v1.v1( LinearIndex ), v.v1.v2( LinearIndex ), v.v2.v1( LinearIndex ), v.v2.v2( LinearIndex ) );
                end
            end
        end

        function v = real( v )
            v.v1.v1 = real( v.v1.v1 );
            v.v1.v2 = real( v.v1.v2 );
            v.v2.v1 = real( v.v2.v1 );
            v.v2.v2 = real( v.v2.v2 );
        end

        function v = imag( v )
            v.v1.v1 = imag( v.v1.v1 );
            v.v1.v2 = imag( v.v1.v2 );
            v.v2.v1 = imag( v.v2.v1 );
            v.v2.v2 = imag( v.v2.v2 );
        end

        function v = conj( v )
            v.v1.v1 = conj( v.v1.v1 );
            v.v1.v2 = conj( v.v1.v2 );
            v.v2.v1 = conj( v.v2.v1 );
            v.v2.v2 = conj( v.v2.v2 );
        end

        function v = angle( v )
            if isreal( v )
                Select = v >= 0;
                v.v1.v1( Select ) = 0;
                v.v1.v2( Select ) = 0;
                v.v2.v1( Select ) = 0;
                v.v2.v2( Select ) = 0;
                v.v1.v1( ~Select ) = QuadDouble.pi.v1.v1;
                v.v1.v2( ~Select ) = QuadDouble.pi.v1.v2;
                v.v2.v1( ~Select ) = QuadDouble.pi.v2.v1;
                v.v2.v2( ~Select ) = QuadDouble.pi.v2.v2;
            else
                v = atan2( imag( v ), real( v ) );
            end
        end

        function [ v, varargout ] = size( v, varargin )
            v = size( v.v1.v1, varargin{:} );
            if nargout > 1
                varargout = num2cell( v( 2:end ) );
                v = v( 1 );
            end
        end

        function v = length( v )
            v = max( size( v ) );
        end

        function v = numel( v )
            v = numel( v.v1.v1 );
        end

        function n = numArgumentsFromSubscript( ~, ~, ~ )
            n = 1;
        end

        function v = end( v, k, n )
            if n == 1
                v = numel( v.v1.v1 );
            else
                v = size( v.v1.v1 );
                if k <= length( v )
                    v = v( k );
                else
                    v = 1;
                end
            end
        end

        function v = repmat( v, varargin )
            v = QuadDouble.Make( repmat( v.v1.v1, varargin{:} ), repmat( v.v1.v2, varargin{:} ), repmat( v.v2.v1, varargin{:} ), repmat( v.v2.v2, varargin{:} ) );
        end

        function v = reshape( v, varargin )
            v.v1.v1 = reshape( v.v1.v1, varargin{:} );
            v.v1.v2 = reshape( v.v1.v2, varargin{:} );
            v.v2.v1 = reshape( v.v2.v1, varargin{:} );
            v.v2.v2 = reshape( v.v2.v2, varargin{:} );
        end

        function v = isequal( a, b, varargin )
            if any( size( a ) ~= size( b ) )
                v = false;
                return
            end
            v = a == b;
            v = all( v( : ) );
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
            v = isempty( v.v1.v1 );
        end

        function v = diag( v, k )
            if nargin < 2
                v = QuadDouble.Make( diag( v.v1.v1 ), diag( v.v1.v2 ), diag( v.v2.v1 ), diag( v.v2.v2 ) );
            else
                v = QuadDouble.Make( diag( v.v1.v1, k ), diag( v.v1.v2, k ), diag( v.v2.v1, k ), diag( v.v2.v2, k ) );
            end
        end

        function v = tril( v, k )
            if nargin < 2
                v = QuadDouble.Make( tril( v.v1.v1 ), tril( v.v1.v2 ), tril( v.v2.v1 ), tril( v.v2.v2 ) );
            else
                v = QuadDouble.Make( tril( v.v1.v1, k ), tril( v.v1.v2, k ), tril( v.v2.v1, k ), tril( v.v2.v2, k ) );
            end
        end

        function v = triu( v, k )
            if nargin < 2
                v = QuadDouble.Make( triu( v.v1.v1 ), triu( v.v1.v2 ), triu( v.v2.v1 ), triu( v.v2.v2 ) );
            else
                v = QuadDouble.Make( triu( v.v1.v1, k ), triu( v.v1.v2, k ), triu( v.v2.v1, k ), triu( v.v2.v2, k ) );
            end
        end

        function v = plus( a, b )
            v = QuadDouble.Plus( a, b );
        end

        function v = minus( a, b )
            v = QuadDouble.Minus( a, b );
        end

        function v = uminus( v )
            v.v1.v1 = -v.v1.v1;
            v.v1.v2 = -v.v1.v2;
            v.v2.v1 = -v.v2.v1;
            v.v2.v2 = -v.v2.v2;
        end

        function v = uplus( v )
        end

        function v = times( a, b )
            v = QuadDouble.Times( a, b );
        end

        function v = mtimes( a, b )
            v = QuadDouble.MTimes( a, b );
        end

        function v = rdivide( a, b )
            v = QuadDouble.RDivide( a, b );
        end

        function v = ldivide( a, b )
            v = QuadDouble.LDivide( a, b );
        end

        function v = mldivide( a, v )
            v = QuadDouble.MLDivide( a, v );
        end

        function v = mrdivide( v, a )
            v = QuadDouble.MRDivide( v, a );
        end

        function v = power( a, b )
            if ~isa( a, 'QuadDouble' )
                a = QuadDouble( a );
            end
            if isa( b, 'QuadDouble' )
                b2 = TimesPowerOf2( b, 2 );
            else
                b2 = 2 * b;
            end
            bHalfInteger = ( double( b2 ) == b2 ) & ( floor( b2 ) == b2 );
            if ~any( bHalfInteger )
                v = exp( b .* log( a ) );
            else
                [ a, b, bHalfInteger ] = QuadDouble.ExpandSingleton( a, b, bHalfInteger );
                v = QuadDouble.ones( size( a ) );
                bNonHalfInteger = find( ~bHalfInteger );
                bHalfInteger = find( bHalfInteger );
                v = QuadDouble.Assign( v, bNonHalfInteger, exp( QuadDouble.Index( b, bNonHalfInteger ) .* log( QuadDouble.Index( a, bNonHalfInteger ) ) ) );
                a = QuadDouble.Index( a, bHalfInteger );
                b = double( QuadDouble.Index( b, bHalfInteger ) );
                Select = find( b < 0 );
                a = QuadDouble.Assign( a, Select, 1 ./ QuadDouble.Index( a, Select ) );
                b = QuadDouble.Assign( b, Select, -b( Select ) );
                vv =  QuadDouble.Index( v, bHalfInteger );
                Select = find( b ~= floor( b ) );
                vv = QuadDouble.Assign( vv, Select, sqrt( QuadDouble.Index( a, Select ) ) );
                Binary = dec2bin( floor( b ) );
                N = size( Binary, 2 );
                Power = a;
                Select = find( Binary( :, end ) == '1' );
                vv = QuadDouble.Assign( vv, Select, QuadDouble.Index( a, Select ) );
                for n = 2 : N
                    Power = Power .* Power;
                    Select = find( Binary( :, end + 1 - n ) == '1' );
                    vv = QuadDouble.Assign( vv, Select, QuadDouble.Index( vv, Select ) .* QuadDouble.Index( Power, Select ) );
                end
                v = QuadDouble.Assign( v, bHalfInteger, vv );
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
                    assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
                    v = v * diag( a .^ d ) / v;
                end
            else
                if nb == 1
                    [ v, d ] = eig( a );
                    d = diag( d );
                    assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
                    v = v * diag( d .^ b ) / v;
                else
                    v = QuadDouble;
                end
            end
        end

        function v = lt( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( ( a2 < b2 ) | ( ( a2 == b2 ) & ( ( a3 < b3 ) | ( ( a3 == b3 ) & ( a4 < b4 ) ) ) ) ) );
        end

        function v = gt( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( ( a2 > b2 ) | ( ( a2 == b2 ) & ( ( a3 > b3 ) | ( ( a3 == b3 ) & ( a4 > b4 ) ) ) ) ) );
        end

        function v = le( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( ( a2 < b2 ) | ( ( a2 == b2 ) & ( ( a3 < b3 ) | ( ( a3 == b3 ) & ( a4 <= b4 ) ) ) ) ) );
        end

        function v = ge( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( ( a2 > b2 ) | ( ( a2 == b2 ) & ( ( a3 > b3 ) | ( ( a3 == b3 ) & ( a4 >= b4 ) ) ) ) ) );
        end

        function v = ne( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 ~= b1 ) | ( a2 ~= b2 ) | ( a3 ~= b3 ) | ( a4 ~= b4 );
        end

        function v = eq( a, b )
            if isa( a, 'QuadDouble' )
                a1 = a.v1.v1; a2 = a.v1.v2; a3 = a.v2.v1; a4 = a.v2.v2;
            else
                a1 = a; a2 = 0; a3 = 0; a4 = 0;
            end
            if isa( b, 'QuadDouble' )
                b1 = b.v1.v1; b2 = b.v1.v2; b3 = b.v2.v1; b4 = b.v2.v2;
            else
                b1 = b; b2 = 0; b3 = 0; b4 = 0;
            end
            v = ( a1 == b1 ) & ( ~isfinite( a1 ) | ( ( a2 == b2 ) & ( a3 == b3 ) & ( a4 == b4 ) ) );
        end

        function v = colon( a, d, b )
            if nargin < 3
                b = d;
                d = 1;
            end
            if ~isa( a, 'QuadDouble' )
                a = QuadDouble( a );
            end
            if ~isa( b, 'QuadDouble' )
                b = QuadDouble( b );
            end
            if ~isa( d, 'QuadDouble' )
                d = QuadDouble( d );
            end
            c = double( floor( ( b - a ) ./ d ) );
            v = a + ( 0:c ) .* d;
        end

        function v = ctranspose( v )
            v.v1.v1 = v.v1.v1';
            v.v1.v2 = v.v1.v2';
            v.v2.v1 = v.v2.v1';
            v.v2.v2 = v.v2.v2';
        end

        function v = transpose( v )
            v.v1.v1 = v.v1.v1.';
            v.v1.v2 = v.v1.v2.';
            v.v2.v1 = v.v2.v1.';
            v.v2.v2 = v.v2.v2.';
        end

        function v = permute( v, DimOrder )
            v.v1.v1 = permute( v.v1.v1, DimOrder );
            v.v1.v2 = permute( v.v1.v2, DimOrder );
            v.v2.v1 = permute( v.v2.v1, DimOrder );
            v.v2.v2 = permute( v.v2.v2, DimOrder );
        end

        function v = ipermute( v, DimOrder )
            v.v1.v1 = ipermute( v.v1.v1, DimOrder );
            v.v1.v2 = ipermute( v.v1.v2, DimOrder );
            v.v2.v1 = ipermute( v.v2.v1, DimOrder );
            v.v2.v2 = ipermute( v.v2.v2, DimOrder );
        end

        function v = horzcat( a, b, varargin )
            if nargin == 1
                v = a;
            elseif nargin > 2
                v = horzcat( horzcat( a, b ), varargin{:} );
            else
                if ~isa( a, 'QuadDouble' )
                    a = QuadDouble( a );
                end
                if ~isa( b, 'QuadDouble' )
                    b = QuadDouble( b );
                end
                x1 = [ a.v1.v1, b.v1.v1 ];
                x2 = [ a.v1.v2, b.v1.v2 ];
                x3 = [ a.v2.v1, b.v2.v1 ];
                x4 = [ a.v2.v2, b.v2.v2 ];
                v = QuadDouble.Make( x1, x2, x3, x4 );
            end
        end

        function v = vertcat( a, b, varargin )
            if nargin == 1
                v = a;
            elseif nargin > 2
                v = vertcat( vertcat( a, b ), varargin{:} );
            else
                if ~isa( a, 'QuadDouble' )
                    a = QuadDouble( a );
                end
                if ~isa( b, 'QuadDouble' )
                    b = QuadDouble( b );
                end
                x1 = [ a.v1.v1; b.v1.v1 ];
                x2 = [ a.v1.v2; b.v1.v2 ];
                x3 = [ a.v2.v1; b.v2.v1 ];
                x4 = [ a.v2.v2; b.v2.v2 ];
                v = QuadDouble.Make( x1, x2, x3, x4 );
            end
        end

        function v = cat( Dim, varargin )
            if Dim == 1
                v = vertcat( varargin{:} );
            elseif Dim == 2
                v = horzcat( varargin{:} );
            else
                % For higher dimensions, convert each to doubles, cat, then rebuild
                x1 = cellfun( @(x) x.v1.v1, varargin, 'UniformOutput', false );
                x2 = cellfun( @(x) x.v1.v2, varargin, 'UniformOutput', false );
                x3 = cellfun( @(x) x.v2.v1, varargin, 'UniformOutput', false );
                x4 = cellfun( @(x) x.v2.v2, varargin, 'UniformOutput', false );
                v = QuadDouble.Make( cat( Dim, x1{:} ), cat( Dim, x2{:} ), cat( Dim, x3{:} ), cat( Dim, x4{:} ) );
            end
        end

        function v = subsref( v, s )
            if strcmp( s(1).type, '.' )
                v = builtin( 'subsref', v, s );
            else
                v.v1.v1 = subsref( v.v1.v1, s );
                v.v1.v2 = subsref( v.v1.v2, s );
                v.v2.v1 = subsref( v.v2.v1, s );
                v.v2.v2 = subsref( v.v2.v2, s );
            end
        end

        function v = subsasgn( v, s, b )
            target_empty = false;
            for k=1:length(s.subs)
                if isempty(s.subs{k})
                    target_empty = true;
                    break;
                end
            end
            if target_empty
                return;
            end

            if ~isa( v, 'QuadDouble' )
                v = QuadDouble( v );
            end
            if ~isa( b, 'QuadDouble' )
                b = QuadDouble( b );
            end
            v.v1.v1 = subsasgn( v.v1.v1, s, b.v1.v1 );
            v.v1.v2 = subsasgn( v.v1.v2, s, b.v1.v2 );
            v.v2.v1 = subsasgn( v.v2.v1, s, b.v2.v1 );
            v.v2.v2 = subsasgn( v.v2.v2, s, b.v2.v2 );
        end

        function v = subsindex( v )
            v = v.v1.v1;
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
            [ v, Indices ] = QuadDouble.Sort( v, Dim, cm, varargin{:} );
        end

        function v = sum( v, Dim )
            if nargin < 2
                Dim = [];
            end
            v = QuadDouble.Sum( v, Dim );
        end

        function v = prod( v, Dim )
            if nargin < 2
                Dim = [];
            end
            v = QuadDouble.Prod( v, Dim );
        end

        function [ v, i ] = max( a, b, Dim )
            if nargin < 3
                Dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if nargout < 2
                v = QuadDouble.Max( a, b, Dim );
            else
                [ v, i ] = QuadDouble.Max( a, b, Dim );
            end
        end

        function [ v, i ] = min( a, b, Dim )
            if nargin < 3
                Dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if nargout < 2
                v = QuadDouble.Min( a, b, Dim );
            else
                [ v, i ] = QuadDouble.Min( a, b, Dim );
            end
        end

        function v = cumsum( v, Dim )
            if nargin < 2
                Dim = [];
            end
            v = QuadDouble.CumSum( v, Dim );
        end

        function v = diff( v, Dim )
            if nargin < 2
                Dim = [];
            end
            v = QuadDouble.Diff( v, Dim );
        end

        function v = cumprod( v, Dim )
            if nargin < 2
                Dim = [];
            end
            v = QuadDouble.CumProd( v, Dim );
        end

        function v = cummax( v, Dim )
            if nargin < 3
                Dim = [];
            end
            v = QuadDouble.CumMax( v, Dim );
        end

        function v = cummin( v, Dim )
            if nargin < 3
                Dim = [];
            end
            v = QuadDouble.CumMin( v, Dim );
        end

        function v = dot( a, b, Dim )
            if nargin < 3
                Dim = [];
            end
            v = QuadDouble.Dot( a, b, Dim );
        end

        function v = norm( v, p )
            if nargin < 2
                v = QuadDouble.Norm( v );
            else
                v = QuadDouble.Norm( v, p );
            end
        end

        function v = abs( v )
            if isreal( v )
                Select = v.v1.v1 < 0;
                v.v1.v1( Select ) = -v.v1.v1( Select );
                v.v1.v2( Select ) = -v.v1.v2( Select );
                v.v2.v1( Select ) = -v.v2.v1( Select );
                v.v2.v2( Select ) = -v.v2.v2( Select );
            else
                real_v = real( v );
                imag_v = imag( v );
                v = sqrt( real_v .* real_v + imag_v .* imag_v );
            end
        end

        function v = sign( v )
            if isreal( v )
                v = sign( v.v1.v1 );
            else
                abs_v = abs( v );
                [ v.v1.v1, v.v1.v2 ] = QuadDouble.QDDivQD( v.v1.v1, v.v1.v2, abs_v.v1.v1, abs_v.v1.v2 );
            end
        end

        function v = floor( v )
            x1 = floor( v.v1.v1 );
            x2 = zeros( size( x1 ) );
            x3 = x2; x4 = x2;
            Select = x1 == v.v1.v1;
            x2( Select ) = floor( v.v1.v2( Select ) );
            Select = Select & ( x2 == v.v1.v2 );
            x3( Select ) = floor( v.v2.v1( Select ) );
            Select = Select & ( x3 == v.v2.v1 );
            x4( Select ) = floor( v.v2.v2( Select ) );
            [ x1, x2, x3, x4 ] = QuadDouble.Renorm4( x1, x2, x3, x4 );
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end
        function v = ceil( v )
            x1 = ceil( v.v1.v1 );
            x2 = zeros( size( x1 ) );
            x3 = x2; x4 = x2;
            Select = x1 == v.v1.v1;
            x2( Select ) = ceil( v.v1.v2( Select ) );
            Select = Select & ( x2 == v.v1.v2 );
            x3( Select ) = ceil( v.v2.v1( Select ) );
            Select = Select & ( x3 == v.v2.v1 );
            x4( Select ) = ceil( v.v2.v2( Select ) );
            [ x1, x2, x3, x4 ] = QuadDouble.Renorm4( x1, x2, x3, x4 );
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end
        function v = fix( v )
            x1 = fix( v.v1.v1 );
            x2 = zeros( size( x1 ) );
            x3 = x2; x4 = x2;
            Select = x1 == v.v1.v1;
            x2( Select ) = fix( v.v1.v2( Select ) );
            Select = Select & ( x2 == v.v1.v2 );
            x3( Select ) = fix( v.v2.v1( Select ) );
            Select = Select & ( x3 == v.v2.v1 );
            x4( Select ) = fix( v.v2.v2( Select ) );
            [ x1, x2, x3, x4 ] = QuadDouble.Renorm4( x1, x2, x3, x4 );
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end
        function v = round( v )
            x1 = round( v.v1.v1 );
            x2 = zeros( size( x1 ) );
            x3 = x2; x4 = x2;
            Select = x1 == v.v1.v1;
            x2( Select ) = round( v.v1.v2( Select ) );
            Select = Select & ( x2 == v.v1.v2 );
            x3( Select ) = round( v.v2.v1( Select ) );
            Select = Select & ( x3 == v.v2.v1 );
            x4( Select ) = round( v.v2.v2( Select ) );
            % This is a simple round that doesn't handle .5 perfectly, but it matches qd_real's nint for most cases.
            [ x1, x2, x3, x4 ] = QuadDouble.Renorm4( x1, x2, x3, x4 );
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end
        function v = realsqrt( v )
            Select = v < 0;
            v.v1.v1( Select ) = NaN;
            v.v1.v2( Select ) = NaN;
            Select = v > 0;
            x = 1 ./ sqrt( v.v1.v1( Select ) );
            vx = v.v1.v1( Select ) .* x;
            t = QuadDouble.Make( v.v1.v1( Select ), v.v1.v2( Select ), v.v2.v1( Select ), v.v2.v2( Select ) ) - QuadDouble.Times( vx, vx );
            t = QuadDouble.Plus( vx, t.v1.v1 .* ( x * 0.5 ) );
            v.v1.v1( Select ) = t.v1.v1;
            v.v1.v2( Select ) = t.v1.v2;
        end

        function v = sqrt( v )
            Select = v ~= 0;
            x = 1 ./ sqrt( v.v1.v1( Select ) );
            vx = v.v1.v1( Select ) .* x;
            t = QuadDouble.Make( v.v1.v1( Select ), v.v1.v2( Select ), v.v2.v1( Select ), v.v2.v2( Select ) ) - QuadDouble.Times( vx, vx );
            t = QuadDouble.Plus( vx, t.v1.v1 .* ( x * 0.5 ) );
            v.v1.v1( Select ) = t.v1.v1;
            v.v1.v2( Select ) = t.v1.v2;
        end

        function v = sqrtm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
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
            % exp( kr + m * log( 2 ) ) = 2^m * exp( r )^k
            % where m and k are integers.  By choosing m appropriately
            % we can make |kr| <= log( 2 ) / 2 = 0.347.  Then exp( r ) is
            % evaluated using the familiar Taylor series.  Reducing the
            % argument substantially speeds up the convergence.
            k = 512.0;
            inv_k = 1.0 / k;
            Threshhold = inv_k .* QuadDouble.eps.v1.v1;

            m = floor( v.v1.v1 ./ QuadDouble.log_2.v1.v1 + 0.5 );
            r = TimesPowerOf2( v - QuadDouble.log_2 .* m, inv_k );

            p = r .* r;
            s = r + TimesPowerOf2( p, 0.5 );
            p = p .* r;
            t = p .* QuadDouble.Make( QuadDouble.InverseFactorial( 3, 1 ), QuadDouble.InverseFactorial( 3, 2 ), QuadDouble.InverseFactorial( 3, 3 ), QuadDouble.InverseFactorial( 3, 4 ) );
            for i = 4 : QuadDouble.NInverseFactorial
                s = s + t;
                p = p .* r;
                t = p .* QuadDouble.Make( QuadDouble.InverseFactorial( i, 1 ), QuadDouble.InverseFactorial( i, 2 ), QuadDouble.InverseFactorial( i, 3 ), QuadDouble.InverseFactorial( i, 4 ) );
                if all( abs( t.v1.v1( : ) ) <= Threshhold )
                    break
                end
            end

            s = s + t;

            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            s = TimesPowerOf2( s, 2.0 ) + s .* s;
            if expm1Flag
                Select = m ~= 0;
                [ s.v1.v1( Select ), s.v1.v2( Select ), s.v2.v1( Select ), s.v2.v2( Select ) ] = QuadDouble.QDPlusDouble( s.v1.v1( Select ), s.v1.v2( Select ), s.v2.v1( Select ), s.v2.v2( Select ), 1.0 );
            else
                s = s + 1.0;
            end

            v = TimesPowerOf2( s, pow2( 1.0, m ) );
            if expm1Flag
                [ v.v1.v1( Select ), v.v1.v2( Select ), v.v2.v1( Select ), v.v2.v2( Select ) ] = QuadDouble.QDPlusDouble( v.v1.v1( Select ), v.v1.v2( Select ), v.v2.v1( Select ), v.v2.v2( Select ), -1.0 );
            end
        end

        function v = expm1( v )
            v = exp( v, true );
        end

        function v = expm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
            v = v * diag( exp( d ) ) / v;
        end

        function x = log( v )
            x = QuadDouble.Make( log( v.v1.v1 ), zeros( size( v.v1.v1 ) ) );
            x = x + v .* exp( -x ) - 1.0;
            x = x + v .* exp( -x ) - 1.0; % slightly paranoid, but does correct e.g. log( exp( QuadDouble( -40 ) ) )
        end

        function v = log2( v )
            v = log( v ) ./ QuadDouble.log_2;
        end

        function v = log10( v )
            v = log( v ) ./ QuadDouble.log_10;
        end

        function v = logm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
            v = v * diag( log( d ) ) / v;
        end

        function v = funm( v, f )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1.v1 ) ) == length( d.v1.v1 ) );
            v = v * diag( f( d ) ) / v;
        end

        function [ sin_v, cos_v ] = sincos( v )
            if ~isreal( v )
                exp_Piv = exp( TimesPowerOf2( v, 1i ) );
                exp_Niv = 1 ./ exp_Piv;
                sin_v = TimesPowerOf2( exp_Piv - exp_Niv, -0.5i );
                cos_v = TimesPowerOf2( exp_Piv + exp_Niv, +0.5  );
                return
            end
            % Strategy.  To compute sin( x ), cos( x ), we choose integers a, b so that
            % x = s + a * ( pi/2 ) + b * ( pi/16 )
            % and |s| <= pi/32.  Using the fact that
            % sin( pi/16 ) = 0.5 * sqrt( 2 - sqrt( 2 + sqrt( 2 ) ) )
            % we can compute sin( x ) from sin( s ), cos( s ).  This greatly increases the convergence of the sine Taylor series.

            z = round( v ./ QuadDouble.piT2 );
            r = v - QuadDouble.piT2 .* z;

            q = floor( r.v1.v1 ./ QuadDouble.piD2.v1.v1 + 0.5 );
            t = r - QuadDouble.piD2 .* q;
            j = q;
            abs_j = abs( j );

            q = floor( t.v1.v1 ./ QuadDouble.piD16.v1.v1 + 0.5 );
            t = t - QuadDouble.piD16 .* q;
            k = q;
            abs_k = abs( k );

            test = ( j >= -2 ) & ( j <= 2 );
            assert( all( test( : ) ) );
            test = abs_k <= 4;
            assert( all( test( : ) ) );

            [ sin_t, cos_t ] = SinCosTaylor( t );

            sin_v = sin_t;
            cos_v = cos_t;

            a = QuadDouble.Make( QuadDouble.CosTable( abs_k + 1, 1 ), QuadDouble.CosTable( abs_k + 1, 2 ), QuadDouble.CosTable( abs_k + 1, 3 ), QuadDouble.CosTable( abs_k + 1, 4 ) );
            b = QuadDouble.Make( QuadDouble.SinTable( abs_k + 1, 1 ), QuadDouble.SinTable( abs_k + 1, 2 ), QuadDouble.SinTable( abs_k + 1, 3 ), QuadDouble.SinTable( abs_k + 1, 4 ) );

            a = reshape( a, size( v ) );
            b = reshape( b, size( v ) );

            a_sin_t = a .* sin_t;
            b_sin_t = b .* sin_t;
            a_cos_t = a .* cos_t;
            b_cos_t = b .* cos_t;

            Select = k > 0;

            [ sin_v.v1.v1( Select ), sin_v.v1.v2( Select ), sin_v.v2.v1( Select ), sin_v.v2.v2( Select ) ] = QuadDouble.QDPlusQD( +a_sin_t.v1.v1( Select ), +a_sin_t.v1.v2( Select ), +a_sin_t.v2.v1( Select ), +a_sin_t.v2.v2( Select ), +b_cos_t.v1.v1( Select ), +b_cos_t.v1.v2( Select ), +b_cos_t.v2.v1( Select ), +b_cos_t.v2.v2( Select ) );
            [ cos_v.v1.v1( Select ), cos_v.v1.v2( Select ), cos_v.v2.v1( Select ), cos_v.v2.v2( Select ) ] = QuadDouble.QDPlusQD( -b_sin_t.v1.v1( Select ), -b_sin_t.v1.v2( Select ), -b_sin_t.v2.v1( Select ), -b_sin_t.v2.v2( Select ), +a_cos_t.v1.v1( Select ), +a_cos_t.v1.v2( Select ), +a_cos_t.v2.v1( Select ), +a_cos_t.v2.v2( Select ) );

            Select = k < 0;

            [ sin_v.v1.v1( Select ), sin_v.v1.v2( Select ), sin_v.v2.v1( Select ), sin_v.v2.v2( Select ) ] = QuadDouble.QDPlusQD( +a_sin_t.v1.v1( Select ), +a_sin_t.v1.v2( Select ), +a_sin_t.v2.v1( Select ), +a_sin_t.v2.v2( Select ), -b_cos_t.v1.v1( Select ), -b_cos_t.v1.v2( Select ), -b_cos_t.v2.v1( Select ), -b_cos_t.v2.v2( Select ) );
            [ cos_v.v1.v1( Select ), cos_v.v1.v2( Select ), cos_v.v2.v1( Select ), cos_v.v2.v2( Select ) ] = QuadDouble.QDPlusQD( +b_sin_t.v1.v1( Select ), +b_sin_t.v1.v2( Select ), +b_sin_t.v2.v1( Select ), +b_sin_t.v2.v2( Select ), +a_cos_t.v1.v1( Select ), +a_cos_t.v1.v2( Select ), +a_cos_t.v2.v1( Select ), +a_cos_t.v2.v2( Select ) );

            Select = j == 1;

            [ sin_v.v1.v1( Select ), sin_v.v1.v2( Select ), sin_v.v2.v1( Select ), sin_v.v2.v2( Select ), cos_v.v1.v1( Select ), cos_v.v1.v2( Select ), cos_v.v2.v1( Select ), cos_v.v2.v2( Select ) ] = deal( +cos_v.v1.v1( Select ), +cos_v.v1.v2( Select ), +cos_v.v2.v1( Select ), +cos_v.v2.v2( Select ), -sin_v.v1.v1( Select ), -sin_v.v1.v2( Select ), -sin_v.v2.v1( Select ), -sin_v.v2.v2( Select ) );

            Select = j == -1;

            [ sin_v.v1.v1( Select ), sin_v.v1.v2( Select ), sin_v.v2.v1( Select ), sin_v.v2.v2( Select ), cos_v.v1.v1( Select ), cos_v.v1.v2( Select ), cos_v.v2.v1( Select ), cos_v.v2.v2( Select ) ] = deal( -cos_v.v1.v1( Select ), -cos_v.v1.v2( Select ), -cos_v.v2.v1( Select ), -cos_v.v2.v2( Select ), +sin_v.v1.v1( Select ), +sin_v.v1.v2( Select ), +sin_v.v2.v1( Select ), +sin_v.v2.v2( Select ) );

            Select = abs_j == 2;

            [ sin_v.v1.v1( Select ), sin_v.v1.v2( Select ), sin_v.v2.v1( Select ), sin_v.v2.v2( Select ), cos_v.v1.v1( Select ), cos_v.v1.v2( Select ), cos_v.v2.v1( Select ), cos_v.v2.v2( Select ) ] = deal( -sin_v.v1.v1( Select ), -sin_v.v1.v2( Select ), -sin_v.v2.v1( Select ), -sin_v.v2.v2( Select ), -cos_v.v1.v1( Select ), -cos_v.v1.v2( Select ), -cos_v.v2.v1( Select ), -cos_v.v2.v2( Select ) );
        end

        function v = sin( v )
            [ v, ~ ] = sincos( v );
        end

        function v = asin( v )
            assert( all( abs( v.v1.v1( : ) ) <= 1 ) );
            v = atan2( v, sqrt( 1 - v.*v ) );
        end

        function v = cos( v )
            [ ~, v ] = sincos( v );
        end

        function v = acos( v )
            assert( all( abs( v.v1.v1( : ) ) <= 1 ) );
            v = atan2( sqrt( 1 - v.*v ), v );
        end

        function v = tan( v )
            [ sin_v, cos_v ] = sincos( v );
            v = sin_v ./ cos_v;
        end

        function v = atan( v )
            v = atan2( v, QuadDouble.Make( 1, 0, 0, 0 ) );
        end

        function v = atan2( y, x )
            r = sqrt( x.*x + y.*y );
            xx = x ./ r;
            yy = y ./ r;
            Select = abs( xx.v1.v1 ) > abs( yy.v1.v1 );
            v = QuadDouble( atan2( y.v1.v1, x.v1.v1 ) );
            [ sin_z, cos_z ] = sincos( v );
            t = yy;
            [ t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ) ] = QuadDouble.QDPlusQD( t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ),      -sin_z.v1.v1( Select ), -sin_z.v1.v2( Select ), -sin_z.v2.v1( Select ), -sin_z.v2.v2( Select ) );
            [ t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ) ] = QuadDouble.QDDivQD( t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ), +cos_z.v1.v1( Select ), +cos_z.v1.v2( Select ), +cos_z.v2.v1( Select ), +cos_z.v2.v2( Select ) );
            Select = ~Select;
            [ t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ) ] = QuadDouble.QDPlusQD( xx.v1.v1( Select ), xx.v1.v2( Select ), xx.v2.v1( Select ), xx.v2.v2( Select ),    -cos_z.v1.v1( Select ), -cos_z.v1.v2( Select ), -cos_z.v2.v1( Select ), -cos_z.v2.v2( Select ) );
            [ t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ) ] = QuadDouble.QDDivQD( t.v1.v1( Select ), t.v1.v2( Select ), t.v2.v1( Select ), t.v2.v2( Select ), -sin_z.v1.v1( Select ), -sin_z.v1.v2( Select ), -sin_z.v2.v1( Select ), -sin_z.v2.v2( Select ) );
            v = v + t;
        end

        function v = sinh( v )
            exp_v = exp( v );
            v = TimesPowerOf2( exp_v - 1 ./ exp_v, 0.5 );
        end

        function v = asinh( v )
            v = log( v + sqrt( v.*v + 1 ) );
        end

        function v = cosh( v )
            exp_v = exp( v );
            v = TimesPowerOf2( exp_v + 1 ./ exp_v, 0.5 );
        end

        function v = acosh( v )
            v = log( v + sqrt( v.*v - 1 ) );
        end

        function [ sinh_v, cosh_v ] = sinhcosh( v )
            exp_Pv = exp( v );
            exp_Nv = 1 ./ exp_Pv;
            sinh_v = TimesPowerOf2( exp_Pv - exp_Nv, 0.5 );
            cosh_v = TimesPowerOf2( exp_Pv + exp_Nv, 0.5 );
        end

        function v = tanh( v )
            [ sinh_v, cosh_v ] = sinhcosh( v );
            v = sinh_v ./ cosh_v;
        end

        function v = atanh( v )
            v = TimesPowerOf2( log( ( 1 + v ) ./ ( 1 - v ) ), 0.5 );
        end

        function v = mod( v, b )
            v = v - b .* floor( v ./ b );
        end

        function v = rem( v, b )
            v = v - b .* fix( v ./ b );
        end

        function [ v, U, p ] = lu( v, type )
            [ m, n ] = size( v );
            p = 1 : m;

            for k = 1 : min( m, n )

                % Find index of largest element below diagonal in k-th column
                [ ~, midx ] = max( abs( QuadDouble.Make( v.v1.v1( k:m, k ), v.v1.v2( k:m, k ), v.v2.v1( k:m, k ), v.v2.v2( k:m, k ) ) ) );
                midx = midx + k - 1;

                % Skip elimination if column is zero
                if v.v1.v1( midx, k ) ~= 0 || v.v1.v2( midx, k ) ~= 0

                    % Swap pivot row
                    if midx ~= k
                        v.v1.v1( [ k midx ], : ) = v.v1.v1( [ midx k ], : );
                        v.v1.v2( [ k midx ], : ) = v.v1.v2( [ midx k ], : );
                        v.v2.v1( [ k midx ], : ) = v.v2.v1( [ midx k ], : );
                        v.v2.v2( [ k midx ], : ) = v.v2.v2( [ midx k ], : );
                        p( [ k midx ] ) = p( [ midx k ] );
                    end

                    % Compute multipliers
                    i = k + 1 : m;
                    [ v.v1.v1( i, k ), v.v1.v2( i, k ), v.v2.v1( i, k ), v.v2.v2( i, k ) ] = QuadDouble.QDDivQD( v.v1.v1( i, k ), v.v1.v2( i, k ), v.v2.v1( i, k ), v.v2.v2( i, k ), v.v1.v1( k, k ), v.v1.v2( k, k ), v.v2.v1( k, k ), v.v2.v2( k, k ) );

                    % Update the remainder of the matrix
                    j = k + 1 : n;
                    % A( i, j ) = A( i, j ) - A( i, k ) .* A( k, j );
                    [ t1, t2, t3, t4 ] = QuadDouble.QDTimesQD( v.v1.v1( i, k ), v.v1.v2( i, k ), v.v2.v1( i, k ), v.v2.v2( i, k ), v.v1.v1( k, j ), v.v1.v2( k, j ), v.v2.v1( k, j ), v.v2.v2( k, j ) );
                    [ v.v1.v1( i, j ), v.v1.v2( i, j ), v.v2.v1( i, j ), v.v2.v2( i, j ) ] = QuadDouble.QDPlusQD( v.v1.v1( i, j ), v.v1.v2( i, j ), v.v2.v1( i, j ), v.v2.v2( i, j ), -t1, -t2, -t3, -t4 );
                end
            end

            if nargout > 1
                % Separate result
                L = tril( v, -1 ) + eye( m, n, 'QuadDouble' );
                U = triu( v );
                if n > m
                    L.v1.v1 = L.v1.v1( :, 1:m );
                    L.v1.v2 = L.v1.v2( :, 1:m );
                    L.v2.v1 = L.v2.v1( :, 1:m );
                    L.v2.v2 = L.v2.v2( :, 1:m );
                elseif n < m
                    U.v1.v1 = U.v1.v1( 1:n, : );
                    U.v1.v2 = U.v1.v2( 1:n, : );
                    U.v2.v1 = U.v2.v1( 1:n, : );
                    U.v2.v2 = U.v2.v2( 1:n, : );
                end
                v = L;

                if nargout > 2
                    if nargin < 2 || ~strcmp( type, 'vector' )
                        pp = eye( m );
                        pp = pp( p, : );
                        p = pp;
                    end
                else
                    invp( p ) = 1 : m;
                    v.v1.v1 = v.v1.v1( invp, : );
                    v.v1.v2 = v.v1.v2( invp, : );
                    v.v2.v1 = v.v2.v1( invp, : );
                    v.v2.v2 = v.v2.v2( invp, : );
                end
            end
        end

        function [ q, v ] = qr( v )
            [ m, n ] = size( v );
            I = eye( m, 'QuadDouble' );
            QT = I;
            for c = 1 : min( m - 1, n )
                x = QuadDouble.Make( v.v1.v1( :, c ), v.v1.v2( :, c ), v.v2.v1( :, c ), v.v2.v2( :, c ) );
                x.v1.v1( 1 : ( c - 1 ) ) = 0;
                x.v1.v2( 1 : ( c - 1 ) ) = 0;
                x.v2.v1( 1 : ( c - 1 ) ) = 0;
                x.v2.v2( 1 : ( c - 1 ) ) = 0;
                alpha = norm( x );
                sign_x_c = sign( QuadDouble.Make( x.v1.v1( c ), x.v1.v2( c ), x.v2.v1( c ), x.v2.v2( c ) ) );
                if sign_x_c ~= 0
                    alpha = -alpha .* sign_x_c;
                end
                a = x;
                [ a.v1.v1( c ), a.v1.v2( c ), a.v2.v1( c ), a.v2.v2( c ) ] = QuadDouble.QDPlusQD( a.v1.v1( c ), a.v1.v2( c ), a.v2.v1( c ), a.v2.v2( c ), -alpha.v1.v1, -alpha.v1.v2, -alpha.v2.v1, -alpha.v2.v2 );
                b = a ./ norm( a );
                if isreal( v )
                    QTn = I - TimesPowerOf2( b .* b.', 2 );
                else
                    QTn = I - ( 1 + ( x' * b ) ./ ( b' * x ) ) .* ( b .* b' );
                end
                QT = QTn * QT;
                v = QTn * v;
            end
            q = QT';
            v = triu( v );
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
                v = QuadDouble.Make( NaN, NaN, NaN, NaN );
            end
        end

        function v = inv( v )
            n = size( v, 1 );
            v = v \ QuadDouble.Make( eye( n ), zeros( n ) );
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
            [ m, n ] = size( v );
            assert( m == n );
            L = QuadDouble.Make( eye( n ), zeros( n ), zeros( n ), zeros( n ) );
            x1 = zeros( 1, n );
            x2 = x1;
            x3 = x1;
            x4 = x1;
            t1 = x1;
            t2 = x1;
            t3 = x1;
            t4 = x1;
            d = QuadDouble.Make( x1, x2, x3, x4 );
            x1( 1 ) = v.v1.v1( 1, 1 );
            x2( 1 ) = v.v1.v2( 1, 1 );
            x3( 1 ) = v.v2.v1( 1, 1 );
            x4( 1 ) = v.v2.v2( 1, 1 );
            d.v1.v1( 1 ) = x1( 1 );
            d.v1.v2( 1 ) = x2( 1 );
            d.v2.v1( 1 ) = x3( 1 );
            d.v2.v2( 1 ) = x4( 1 );
            idxs = 2 : n;
            [ L.v1.v1( idxs, 1 ), L.v1.v2( idxs, 1 ), L.v2.v1( idxs, 1 ), L.v2.v2( idxs, 1 ) ] = QuadDouble.QDDivQD( v.v1.v1( idxs, 1 ), v.v1.v2( idxs, 1 ), v.v2.v1( idxs, 1 ), v.v2.v2( idxs, 1 ), x1( 1 ), x2( 1 ), x3( 1 ), x4( 1 ) );
            for j = 2 : n
                idxs = 1 : j - 1;
                [ x1( idxs ), x2( idxs ), x3( idxs ), x4( idxs ) ] = QuadDouble.QDTimesQD( conj( L.v1.v1( j, idxs ) ), conj( L.v1.v2( j, idxs ) ), conj( L.v2.v1( j, idxs ) ), conj( L.v2.v2( j, idxs ) ), d.v1.v1( idxs ), d.v1.v2( idxs ), d.v2.v1( idxs ), d.v2.v2( idxs ) );
                [ t1( idxs ), t2( idxs ), t3( idxs ), t4( idxs ) ] = QuadDouble.QDTimesQD( L.v1.v1( j, idxs ), L.v1.v2( j, idxs ), L.v2.v1( j, idxs ), L.v2.v2( j, idxs ), x1( idxs ), x2( idxs ), x3( idxs ), x4( idxs ) );
                t = sum( QuadDouble.Make( t1( idxs ), t2( idxs ), t3( idxs ), t4( idxs ) ) );
                [ x1( j ), x2( j ), x3( j ), x4( j ) ] = QuadDouble.QDPlusQD( v.v1.v1( j, j ), v.v1.v2( j, j ), v.v2.v1( j, j ), v.v2.v2( j, j ), -t.v1.v1, -t.v1.v2, -t.v2.v1, -t.v2.v2 );
                d.v1.v1( j ) = x1( j );
                d.v1.v2( j ) = x2( j );
                d.v2.v1( j ) = x3( j );
                d.v2.v2( j ) = x4( j );
                if j < n
                    jdxs = j + 1 : n;
                    [ s1, s2, s3, s4 ] = QuadDouble.QDTimesQD( L.v1.v1( jdxs, idxs ), L.v1.v2( jdxs, idxs ), L.v2.v1( jdxs, idxs ), L.v2.v2( jdxs, idxs ), x1( idxs ), x2( idxs ), x3( idxs ), x4( idxs ) );
                    tt = sum( QuadDouble.Make( s1, s2, s3, s4 ), 2 );
                    [ t1( jdxs ), t2( jdxs ), t3( jdxs ), t4( jdxs ) ] = QuadDouble.QDPlusQD( v.v1.v1( jdxs, j ), v.v1.v2( jdxs, j ), v.v2.v1( jdxs, j ), v.v2.v2( jdxs, j ), -tt.v1.v1, -tt.v1.v2, -tt.v2.v1, -tt.v2.v2 );
                    [ L.v1.v1( jdxs, j ), L.v1.v2( jdxs, j ), L.v2.v1( jdxs, j ), L.v2.v2( jdxs, j ) ] = QuadDouble.QDDivQD( t1( jdxs ), t2( jdxs ), t3( jdxs ), t4( jdxs ), x1( j ), x2( j ), x3( j ), x4( j ) );
                end
            end
            if nargin < 2 || ~strcmp( type, 'vector_d' )
                d = diag( d );
            else
                d = d.';
            end
        end

        function [ v, d ] = eig( x )
            [ v, d ] = eig( x.v1.v1 );
            v = QuadDouble( v );
            d = QuadDouble( diag( d ) );
            C = length( d );
            I = eye( C, 'QuadDouble' );
            for c = 1 : C
                vi = QuadDouble.Make( v.v1.v1( :, c ), v.v1.v2( :, c ), v.v2.v1( :, c ), v.v2.v2( :, c ) );
                dii = QuadDouble.Make( d.v1.v1( c, 1 ), d.v1.v2( c, 1 ), d.v2.v1( c, 1 ), d.v2.v2( c, 1 ) );
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
                        % err = oerr;
                        dii = odii;
                        break
                    end
                    if ( err == 0 ) || ( err == oerr )
                        break
                    end
                end
                % disp( err );
                d.v1.v1( c, 1 ) = dii.v1.v1;
                d.v1.v2( c, 1 ) = dii.v1.v2;
                d.v2.v1( c, 1 ) = dii.v2.v1;
                d.v2.v2( c, 1 ) = dii.v2.v2;
                v.v1.v1( :, c ) = vi.v1.v1;
                v.v1.v2( :, c ) = vi.v1.v2;
                v.v2.v1( :, c ) = vi.v2.v1;
                v.v2.v2( :, c ) = vi.v2.v2;
            end
            if nargout < 2
                v = d;
            else
                d = diag( d );
            end
        end

        function w = conv( u, v )

            RowVector = size( u, 1 ) == 1 && size( v, 1 ) == 1;

            u = Vec( u );
            v = Vec( v );

            M = size( u, 1 );
            N = size( v, 1 );

            K = M + N - 1;

            w = QuadDouble.zeros( K, 1 );

            for k = 1 : K

                j = max( 1, k + 1 - N ) : min( k, M );
                i = k - j + 1;

                wk = QuadDouble.Dot( QuadDouble.Make( u.v1.v1( j ), u.v1.v2( j ), u.v2.v1( j ), u.v2.v2( j ) ), QuadDouble.Make( v.v1.v1( i ), v.v1.v2( i ), v.v2.v1( i ), v.v2.v2( i ) ) );
                w.v1.v1( k ) = wk.v1.v1;
                w.v1.v2( k ) = wk.v1.v2;
                w.v2.v1( k ) = wk.v2.v1;
                w.v2.v2( k ) = wk.v2.v2;

            end

            if RowVector
                w = w.';
            end

        end

        function [ C, ia, ic ] = unique( A, varargin )
            Rows = strcmpi( varargin, 'rows' );
            if any( Rows )
                varargin( Rows ) = [];
                RowFlag = 0;
            else
                RowFlag = size( A, 1 ) == 1;
                A = Vec( A );
            end
            Size = size( A );
            A = [ A.v1.v1, A.v1.v2, A.v2.v1, A.v2.v2 ];
            [ C, ia, ic ] = unique( A, 'rows', varargin{:} );
            n = Size( 2 );
            C = QuadDouble.Make( C( :, 1 : n ), C( :, ( n + 1 ) : ( 2 * n ) ), C( :, ( 2 * n ) + 1 : ( 3 * n ) ), C( :, ( 3 * n ) + 1 : ( 4 * n ) ) );
            if RowFlag
                C = C.';
            end
        end

        function v = mean( v, Dim )
            if ( nargin < 2 ) || isempty( Dim )
                Dim = find( size( v.v1.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            elseif strcmpi( Dim, 'all' )
                v = Vec( v );
                Dim = 1;
            end
            Size = size( v.v1.v1 );
            n = prod( Size( Dim ) );
            v = sum( v, Dim ) ./ n;
        end

        function v = median( v, Dim )
            if ( nargin < 2 ) || isempty( Dim )
                Dim = find( size( v.v1.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            elseif strcmpi( Dim, 'all' )
                v = Vec( v );
                Dim = 1;
            end
            Size = size( v.v1.v1 );
            n = prod( Size( Dim ) );

            if n == 0
                Size( Dim ) = 0;
                v = QuadDouble.Make( NaN( Size ), NaN( Size ) );
                return;
            end

            NotDim = setdiff( 1 : numel( Size ), Dim );
            v = reshape( permute( v, [ Dim, NotDim ] ), [ n, Size( NotDim ) ] );

            v = sort( v, 1 );

            if mod( n, 2 ) == 1
                Middle = ( n + 1 ) * 0.5;
                v = QuadDouble.Make( v.v1.v1( Middle, : ), v.v1.v2( Middle, : ), v.v2.v1( Middle, : ), v.v2.v2( Middle, : ) );
            else
                Middle1 = n * 0.5;
                Middle2 = Middle1 + 1;
                m1 = QuadDouble.Make( v.v1.v1( Middle1, : ), v.v1.v2( Middle1, : ), v.v2.v1( Middle1, : ), v.v2.v2( Middle1, : ) );
                m2 = QuadDouble.Make( v.v1.v1( Middle2, : ), v.v1.v2( Middle2, : ), v.v2.v1( Middle2, : ), v.v2.v2( Middle2, : ) );
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
                Dim = find( size( v.v1.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            end

            Size = size( v.v1.v1 );
            n = prod( Size( Dim ) );

            if ( n == 0 ) || ( ( n == 1 ) && ( Flag == 0 ) )
                Size( Dim ) = 1;
                v = QuadDouble.Make( NaN( Size ), NaN( Size ) );
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
            [ X1, Y1 ] = meshgrid( x.v1.v1, y.v1.v1 );
            [ X2, Y2 ] = meshgrid( x.v1.v2, y.v1.v2 );
            [ X3, Y3 ] = meshgrid( x.v2.v1, y.v2.v1 );
            [ X4, Y4 ] = meshgrid( x.v2.v2, y.v2.v2 );
            X = QuadDouble.Make( X1, X2, X3, X4 );
            Y = QuadDouble.Make( Y1, Y2, Y3, Y4 );
        end

        function y = linspace( a, b, n )
            if nargin < 3
                n = 100;
            end

            if ~isa( a, 'QuadDouble' )
                a = QuadDouble( a );
            end

            if ~isa( b, 'QuadDouble' )
                b = QuadDouble( b );
            end

            if n < 1
                y = QuadDouble.zeros( 0, 1 );
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
                y.v1.v1( end ) = b.v1.v1;
                y.v1.v2( end ) = b.v1.v2;
                y.v2.v1( end ) = b.v2.v1;
                y.v2.v2( end ) = b.v2.v2;
            end
        end

    end

    methods ( Static )
        function v = IsEqualWithExpansion( a, b, varargin )
            v = a == b;
            v = all( v( : ) );
            if nargin > 2
                for i = 1 : length( varargin )
                    if ~v
                        break
                    end
                    v = v && IsEqualWithExpansion( a, varargin{i} );
                end
            end
        end

        function v = ones( varargin )
            v = QuadDouble.Make( ones( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = zeros( varargin )
            v = QuadDouble.Make( zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = eye( varargin )
            v = QuadDouble.Make( eye( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = nan( varargin )
            v = QuadDouble.Make( nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ) );
        end

        function v = inf( varargin )
            v = QuadDouble.Make( inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ) );
        end

        function v = rand( varargin )
            t = rand( varargin{:}, 'double' );
            v = QuadDouble.Make( t, eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ), eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ), eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ) );
        end

        function v = randn( varargin )
            t = randn( varargin{:}, 'double' );
            v = QuadDouble.Make( t, eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ), eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ), eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ) );
        end

        function v = randi( imax, varargin )
            v = QuadDouble.Make( randi( imax, varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = Plus( a, b )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            if isa( a, 'QuadDouble' )
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusDouble( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, double( b ) );
                end
            else
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusDouble( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, double( a ) );
                else
                    [ x1, x2 ] = QuadDouble.TwoSum( double( a ), double( b ) ); x3 = zeros( size( x1 ) ); x4 = zeros( size( x1 ) );
                end
            end
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end

        function v = Minus( a, b )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            if isa( a, 'QuadDouble' )
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, -b.v1.v1, -b.v1.v2, -b.v2.v1, -b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusDouble( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, -double( b ) );
                end
            else
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDPlusDouble( -b.v1.v1, -b.v1.v2, -b.v2.v1, -b.v2.v2, double( a ) );
                else
                    [ x1, x2 ] = QuadDouble.TwoSum( double( a ), -double( b ) ); x3 = zeros( size( x1 ) ); x4 = zeros( size( x1 ) );
                end
            end
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end

        function v = Times( a, b )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            if isa( a, 'QuadDouble' )
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.QDTimesDouble( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, double( b ) );
                end
            else
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDTimesDouble( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, double( a ) );
                else
                    [ x1, x2 ] = QuadDouble.TwoProd( double( a ), double( b ) ); x3 = zeros( size( x1 ) ); x4 = zeros( size( x1 ) );
                end
            end
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end

        function v = MTimes( a, b )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            [ R, c ] = size( a );
            [ r, C ] = size( b );
            if ( ( R == 1 ) && ( c == 1 ) ) || ( ( r == 1 ) && ( C == 1 ) )
                v = QuadDouble.Times( a, b );
                return
            end
            v = QuadDouble.Make( zeros( R, C ), zeros( R, C ), zeros( R, C ), zeros( R, C ) );
            if isa( b, 'QuadDouble' )
                for c = 1 : C
                    t = QuadDouble.Sum( a .* QuadDouble.Make( b.v1.v1( :, c ).', b.v1.v2( :, c ).', b.v2.v1( :, c ).', b.v2.v2( :, c ).' ), 2 );
                    v.v1.v1( :, c ) = t.v1.v1;
                    v.v1.v2( :, c ) = t.v1.v2;
                    v.v2.v1( :, c ) = t.v2.v1;
                    v.v2.v2( :, c ) = t.v2.v2;
                end
            else
                for c = 1 : C
                    t = QuadDouble.Sum( a .* b( :, c ).', 2 );
                    v.v1.v1( :, c ) = t.v1.v1;
                    v.v1.v2( :, c ) = t.v1.v2;
                    v.v2.v1( :, c ) = t.v2.v1;
                    v.v2.v2( :, c ) = t.v2.v2;
                end
            end
        end

        function v = RDivide( a, b )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            if isa( a, 'QuadDouble' )
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivDouble( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, double( b ) );
                end
            else
                if isa( b, 'QuadDouble' )
                    da = double( a );
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivQD( da, zeros( size( da ) ), zeros( size( da ) ), zeros( size( da ) ), b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.DoubleDivDouble( double( a ), double( b ) );
                end
            end
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end

        function v = LDivide( b, a )
            if isa( a, 'BaseExtDouble' ) && ~isa( a, 'QuadDouble' ); a = QuadDouble( a ); end
            if isa( b, 'BaseExtDouble' ) && ~isa( b, 'QuadDouble' ); b = QuadDouble( b ); end
            if isa( a, 'QuadDouble' )
                if isa( b, 'QuadDouble' )
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivDouble( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, double( b ) );
                end
            else
                if isa( b, 'QuadDouble' )
                    da = double( a );
                    [ x1, x2, x3, x4 ] = QuadDouble.QDDivQD( da, zeros( size( da ) ), zeros( size( da ) ), zeros( size( da ) ), b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                else
                    [ x1, x2, x3, x4 ] = QuadDouble.DoubleDivDouble( double( a ), double( b ) );
                end
            end
            v = QuadDouble.Make( x1, x2, x3, x4 );
        end

        function v = MLDivide( a, v )
            [ Ra, Ca ] = size( a );
            [ Rv, Cv ] = size( v );
            if ( ( Ra == 1 ) && ( Ca == 1 ) ) || ( ( Rv == 1 ) && ( Cv == 1 ) )
                v = QuadDouble.LDivide( a, v );
                return
            end
            if ~isa( v, 'QuadDouble' )
                v = QuadDouble( v );
            end
            assert( Ra == Rv );
            if Ra ~= Ca
                [ q, a ] = qr( a );
                v = q' * v;
                v = BackSubstitution( v, a );
                return
            end
            if QuadDouble.IsEqualWithExpansion( triu( a, 1 ), 0 )
                % Lower triangular
                v = ForwardElimination( v, a );
                return
            elseif QuadDouble.IsEqualWithExpansion( tril( a, -1 ), 0 )
                % Upper triangular
                v = BackSubstitution( v, a );
                return
            elseif QuadDouble.IsEqualWithExpansion( a, a' )
                [ L, d ] = ldl( a, 'vector_d' );
                if all( all( isfinite( L ) ) ) && all( isfinite( d ) )
                    % Positive definite
                    v = ForwardElimination( v, L );
                    v = v ./ d;
                    v = BackSubstitution( v, L' );
                    return
                end
            end
            % Triangular factorization
            [ L, U, p ] = lu( a, 'vector' );

            % Permutation and forward elimination
            v.v1.v1 = v.v1.v1( p, : );
            v.v1.v2 = v.v1.v2( p, : );
            v.v2.v1 = v.v2.v1( p, : );
            v.v2.v2 = v.v2.v2( p, : );
            v = ForwardElimination( v, L );

            % Back substitution
            v = BackSubstitution( v, U );
        end

        function v = MRDivide( v, a )
            v = QuadDouble.MLDivide( a', v' )';
        end

        function [ v, Indices ] = Sort( v, Dim, cm, varargin )
            if nargin < 2 || isempty( Dim )
                Dim = find( size( v.v1.v1 ) > 1, 1 );
                if isempty( Dim )
                    Dim = 1;
                end
            end
            if nargin < 3 || isempty( cm )
                cm = 'auto';
            end
            if isa( v, 'QuadDouble' )
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
                Size = size( v.v1.v1 );
                if any( Size == 0 )
                    Indices = [];
                    return
                end
                Blocks = arrayfun( @( x ) ones( x, 1 ), Size, 'UniformOutput', false );
                Blocks{ Dim } = Size( Dim );
                xv1 = mat2cell( v.v1.v1, Blocks{:} );
                xv2 = mat2cell( v.v1.v2, Blocks{:} );
                xa1 = mat2cell( a.v1.v1, Blocks{:} );
                xa2 = mat2cell( a.v1.v2, Blocks{:} );
                xb1 = mat2cell( b.v1.v1, Blocks{:} );
                xb2 = mat2cell( b.v1.v2, Blocks{:} );
                Indices = cell( size( xv1 ) );
                for i = 1 : numel( xv1 )
                    [ ~, Indices{ i } ] = sortrows( [ xa1{ i }( : ), xa2{ i }( : ), xb1{ i }( : ), xb2{ i }( : ) ], varargin{:} );
                    xv1{ i } = xv1{ i }( Indices{ i } );
                    xv2{ i } = xv2{ i }( Indices{ i } );
                end
                Indices = cell2mat( Indices );
                v       = QuadDouble.Make( cell2mat( xv1 ), cell2mat( xv2 ) );
            else
                if nargout > 1
                    [ v, Indices ] = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                else
                    v = sort( v, Dim, 'ComparisonMethod', cm, varargin{:} );
                end
                v = QuadDouble( v );
            end
        end

        function s = Sum( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                for i = 2 : Length
                    s = QuadDouble.Plus( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
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
                    s = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = QuadDouble.Plus( s, x{ i } );
                end
            end
        end

        function c = CumSum( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1.v1;
                c2{1} = s.v1.v2;
                for i = 2 : Length
                    s = QuadDouble.Plus( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
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
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
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
                    s = QuadDouble.Plus( s, x{ i } );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
                end
            end
            c = QuadDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function c = Diff( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } );
                    d = QuadDouble.Minus( t, s );
                    c1{i} = d.v1.v1;
                    c2{i} = d.v1.v2;
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
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
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
                    d = QuadDouble.Minus( t, s );
                    c1{i} = d.v1.v1;
                    c2{i} = d.v1.v2;
                    s = t;
                end
            end
            c = QuadDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function s = Prod( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = QuadDouble.Make( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                for i = 2 : Length
                    s = QuadDouble.Times( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
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
                    s = QuadDouble.Make( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = QuadDouble.Times( s, x{ i } );
                end
            end
        end

        function c = CumProd( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1.v1;
                c2{1} = s.v1.v2;
                for i = 2 : Length
                    s = QuadDouble.Times( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
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
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
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
                    s = QuadDouble.Times( s, x{ i } );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
                end
            end
            c = QuadDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
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
                    v = QuadDouble;
                end
            end
        end

        function [ s, i ] = Max( a, b, Dim )
            if isempty( b )
                if isempty( a )
                    s = QuadDouble;
                    i = [];
                    return
                end
                if isa( a, 'QuadDouble' )
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a.v1.v1 ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a.v1.v1 );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1.v1, Blocks{:} );
                    x2 = mat2cell( a.v1.v2, Blocks{:} );
                    x3 = mat2cell( a.v2.v1, Blocks{:} );
                    x4 = mat2cell( a.v2.v2, Blocks{:} );
                    s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = QuadDouble.Max( QuadDouble.Make( x1{ j }, x2{ j }, x3{ j }, x4{ j } ), s );
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
                        s = QuadDouble.Max( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'QuadDouble' )
                    a = QuadDouble( a );
                end
                if ~isa( b, 'QuadDouble' )
                    b = QuadDouble( b );
                end
                [ a, b ] = QuadDouble.ExpandSingleton( a, b );
                i = ( a.v1.v1 > b.v1.v1 ) | ( ( a.v1.v1 == b.v1.v1 ) & ( a.v1.v2 > b.v1.v2 ) );
                s = b;
                s.v1.v1( i ) = a.v1.v1( i );
                s.v1.v2( i ) = a.v1.v2( i );
                s.v2.v1( i ) = a.v2.v1( i );
                s.v2.v2( i ) = a.v2.v2( i );
            end
        end

        function c = CumMax( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1.v1;
                c2{1} = s.v1.v2;
                for i = 2 : Length
                    s = QuadDouble.Max( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
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
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
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
                    s = QuadDouble.Max( s, x{ i } );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
                end
            end
            c = QuadDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function [ s, i ] = Min( a, b, Dim )
            if isempty( b )
                if isa( a, 'QuadDouble' )
                    if nargin < 3 || isempty( Dim )
                        Dim = find( size( a.v1.v1 ) > 1, 1 );
                        if isempty( Dim )
                            Dim = 1;
                        end
                    end
                    Size = size( a.v1.v1 );
                    Length = Size( Dim );
                    Blocks = num2cell( Size );
                    Blocks{ Dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1.v1, Blocks{:} );
                    x2 = mat2cell( a.v1.v2, Blocks{:} );
                    x3 = mat2cell( a.v2.v1, Blocks{:} );
                    x4 = mat2cell( a.v2.v2, Blocks{:} );
                    s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = QuadDouble.Min( QuadDouble.Make( x1{ j }, x2{ j }, x3{ j }, x4{ j } ), s );
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
                        s = QuadDouble.Min( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'QuadDouble' )
                    a = QuadDouble( a );
                end
                if ~isa( b, 'QuadDouble' )
                    b = QuadDouble( b );
                end
                [ a, b ] = QuadDouble.ExpandSingleton( a, b );
                i = ( a.v1.v1 < b.v1.v1 ) | ( ( a.v1.v1 == b.v1.v1 ) & ( a.v1.v2 < b.v1.v2 ) );
                s = b;
                s.v1.v1( i ) = a.v1.v1( i );
                s.v1.v2( i ) = a.v1.v2( i );
                s.v2.v1( i ) = a.v2.v1( i );
                s.v2.v2( i ) = a.v2.v2( i );
            end
        end

        function c = CumMin( v, Dim )
            if isa( v, 'QuadDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1.v1, Blocks{:} );
                x2 = mat2cell( v.v1.v2, Blocks{:} );
                x3 = mat2cell( v.v2.v1, Blocks{:} );
                x4 = mat2cell( v.v2.v2, Blocks{:} );
                s = QuadDouble.Make( x1{ 1 }, x2{ 1 }, x3{ 1 }, x4{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1.v1;
                c2{1} = s.v1.v2;
                for i = 2 : Length
                    s = QuadDouble.Min( s, QuadDouble.Make( x1{ i }, x2{ i }, x3{ i }, x4{ i } ) );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
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
                    c = QuadDouble.Make( zeros( Size ), zeros( Size ) );
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
                    s = QuadDouble.Min( s, x{ i } );
                    c1{i} = s.v1.v1;
                    c2{i} = s.v1.v2;
                end
            end
            c = QuadDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = Dot( a, b, Dim )
            if nargin < 3
                Dim = [];
            end
            if ( length( a ) == numel( a ) ) && ( length( b ) == numel( b ) )
                a = Vec( a );
                b = Vec( b );
            end
            v = QuadDouble.Sum( QuadDouble.Times( a, b ), Dim );
        end

        function [ varargout ] = ExpandSingleton( varargin )
            l = cellfun( @( x ) length( size( x ) ), varargin, 'UniformOutput', true );
            n = length( varargin );
            ss = ones( n, max( l ) );
            for i = 1 : n
                ss( i, 1 : l( i ) ) = size( varargin{ i } );
            end
            s = max( ss, [], 1 );
            varargout = cell( 1, n );
            for i = 1 : n
                varargout{ i } = repmat( varargin{ i }, s ./ ss( i, : ) );
            end
        end
    end

    methods ( Access = private )
        function v = Vec( v )
            v.v1.v1 = v.v1.v1( : );
            v.v1.v2 = v.v1.v2( : );
            v.v2.v1 = v.v2.v1( : );
            v.v2.v2 = v.v2.v2( : );
        end

        function v = TimesPowerOf2( v, b )
            assert( isa( b, 'double' ) );
            v.v1.v1 = v.v1.v1 .* b;
            v.v1.v2 = v.v1.v2 .* b;
            v.v2.v1 = v.v2.v1 .* b;
            v.v2.v2 = v.v2.v2 .* b;
        end

        function v = ForwardElimination( v, L )
            % For lower triangular L, x = ForwardElimination( b, L ) solves L*x = b.
            [ m, n ] = size( L );
            mn = min( m, n );
            [ vm, vn ] = size( v );
            if vm < n
                v = [ v; zeros( n - vm, vn ) ];
            elseif vm > n
                v.v1.v1( ( n+1 ):vm, : ) = [];
                v.v1.v2( ( n+1 ):vm, : ) = [];
                v.v2.v1( ( n+1 ):vm, : ) = [];
                v.v2.v2( ( n+1 ):vm, : ) = [];
            end
            if isa( L, 'QuadDouble' )
                [ v.v1.v1( 1, : ), v.v1.v2( 1, : ), v.v2.v1( 1, : ), v.v2.v2( 1, : ) ] = QuadDouble.QDDivQD( v.v1.v1( 1, : ), v.v1.v2( 1, : ), v.v2.v1( 1, : ), v.v2.v2( 1, : ), L.v1.v1( 1, 1 ), L.v1.v2( 1, 1 ), L.v2.v1( 1, 1 ), L.v2.v2( 1, 1 ) );
                for k = 2 : mn
                    j = 1 : k - 1;
                    [ t1, t2, t3, t4 ] = QuadDouble.QDTimesQD( v.v1.v1( j, : ), v.v1.v2( j, : ), v.v2.v1( j, : ), v.v2.v2( j, : ), L.v1.v1( k, j ).', L.v1.v2( k, j ).', L.v2.v1( k, j ).', L.v2.v2( k, j ).' );
                    t = QuadDouble.Sum( QuadDouble.Make( t1, t2, t3, t4 ), 1 );
                    [ t1, t2, t3, t4 ] = QuadDouble.QDPlusQD( v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ), -t.v1.v1, -t.v1.v2, -t.v2.v1, -t.v2.v2 );
                    [ v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ) ] = QuadDouble.QDDivQD( t1, t2, t3, t4, L.v1.v1( k, k ), L.v1.v2( k, k ), L.v2.v1( k, k ), L.v2.v2( k, k ) );
                end
            else
                [ v.v1.v1( 1, : ), v.v1.v2( 1, : ), v.v2.v1( 1, : ), v.v2.v2( 1, : ) ] = QuadDouble.QDDivDouble( v.v1.v1( 1, : ), v.v1.v2( 1, : ), v.v2.v1( 1, : ), v.v2.v2( 1, : ), L( 1, 1 ) );
                for k = 2 : mn
                    j = 1 : k - 1;
                    [ t1, t2, t3, t4 ] = QuadDouble.QDTimesDouble( v.v1.v1( j, : ), v.v1.v2( j, : ), v.v2.v1( j, : ), v.v2.v2( j, : ), L( k, j ).' );
                    t = QuadDouble.Sum( QuadDouble.Make( t1, t2, t3, t4 ), 1 );
                    [ t1, t2, t3, t4 ] = QuadDouble.QDPlusQD( v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ), -t.v1.v1, -t.v1.v2, -t.v2.v1, -t.v2.v2 );
                    [ v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ) ] = QuadDouble.QDDivDouble( t1, t2, t3, t4, L( k, k ) );
                end
            end
        end

        function v = BackSubstitution( v, U )
            % For upper triangular U, x = BackSubstitution( b, U ) solves U*x = b.
            [ m, n ] = size( U );
            mn = min( m, n );
            [ vm, vn ] = size( v );
            if vm < n
                v = [ v; zeros( n - vm, vn ) ];
            elseif vm > n
                v.v1.v1( ( n+1 ):vm, : ) = [];
                v.v1.v2( ( n+1 ):vm, : ) = [];
                v.v2.v1( ( n+1 ):vm, : ) = [];
                v.v2.v2( ( n+1 ):vm, : ) = [];
            end
            if isa( U, 'QuadDouble' )
                [ v.v1.v1( mn, : ), v.v1.v2( mn, : ), v.v2.v1( mn, : ), v.v2.v2( mn, : ) ] = QuadDouble.QDDivQD( v.v1.v1( mn, : ), v.v1.v2( mn, : ), v.v2.v1( mn, : ), v.v2.v2( mn, : ), U.v1.v1( mn, mn ), U.v1.v2( mn, mn ), U.v2.v1( mn, mn ), U.v2.v2( mn, mn ) );
                for k = mn - 1 : -1 : 1
                    j = k + 1 : n;
                    [ t1, t2, t3, t4 ] = QuadDouble.QDTimesQD( v.v1.v1( j, : ), v.v1.v2( j, : ), v.v2.v1( j, : ), v.v2.v2( j, : ), U.v1.v1( k, j ).', U.v1.v2( k, j ).', U.v2.v1( k, j ).', U.v2.v2( k, j ).' );
                    t = QuadDouble.Sum( QuadDouble.Make( t1, t2, t3, t4 ), 1 );
                    [ t1, t2, t3, t4 ] = QuadDouble.QDPlusQD( v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ), -t.v1.v1, -t.v1.v2, -t.v2.v1, -t.v2.v2 );
                    [ v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ) ] = QuadDouble.QDDivQD( t1, t2, t3, t4, U.v1.v1( k, k ), U.v1.v2( k, k ), U.v2.v1( k, k ), U.v2.v2( k, k ) );
                end
            else
                [ v.v1.v1( mn, : ), v.v1.v2( mn, : ), v.v2.v1( mn, : ), v.v2.v2( mn, : ) ] = QuadDouble.QDDivDouble( v.v1.v1( mn, : ), v.v1.v2( mn, : ), v.v2.v1( mn, : ), v.v2.v2( mn, : ), U( mn, mn ) );
                for k = mn - 1 : -1 : 1
                    j = k + 1 : n;
                    [ t1, t2, t3, t4 ] = QuadDouble.QDTimesDouble( v.v1.v1( j, : ), v.v1.v2( j, : ), v.v2.v1( j, : ), v.v2.v2( j, : ), U( k, j ).' );
                    t = QuadDouble.Sum( QuadDouble.Make( t1, t2, t3, t4 ), 1 );
                    [ t1, t2, t3, t4 ] = QuadDouble.QDPlusQD( v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ), -t.v1.v1, -t.v1.v2, -t.v2.v1, -t.v2.v2 );
                    [ v.v1.v1( k, : ), v.v1.v2( k, : ), v.v2.v1( k, : ), v.v2.v2( k, : ) ] = QuadDouble.QDDivDouble( t1, t2, t3, t4, U( k, k ) );
                end
            end
        end

        function v = SinTaylor( v )
            Threshhold = 0.5 .* abs( v.v1.v1( : ) ) .* QuadDouble.eps.v1.v1;
            x = - v .* v;
            r = v;
            for i = 3 : 2 : QuadDouble.NInverseFactorial
                r = r .* x;
                t = r .* QuadDouble.Make( QuadDouble.InverseFactorial( i, 1 ), QuadDouble.InverseFactorial( i, 2 ), QuadDouble.InverseFactorial( i, 3 ), QuadDouble.InverseFactorial( i, 4 ) );
                v = v + t;
                if all( abs( t.v1.v1( : ) ) <= Threshhold )
                    break
                end
            end
        end

        function [ sin_v, cos_v ] = SinCosTaylor( v )
            sin_v = SinTaylor( v );
            cos_v = sqrt( 1 - sin_v .* sin_v );
        end
    end

    methods ( Static, Access = { ?QuadDouble, ?OctDouble } )
        function v = Make( a1, a2, a3, a4 )
            v = QuadDouble;
            v.v1.v1 = a1;
            if nargin < 2
                a2 = zeros( size( a1 ) );
            end
            v.v1.v2 = a2;
            if nargin < 3
                a3 = zeros( size( a1 ) );
            end
            v.v2.v1 = a3;
            if nargin < 4
                a4 = zeros( size( a1 ) );
            end
            v.v2.v2 = a4;
        end

        function [ s1, s2 ] = Normalize( a1, a2 )
            s1 = a1 + a2;
            t = s1 - a1;
            s2 = a2 - t;
        end



        function [ s1, s2 ] = QuickTwoSum( a, b )
            s1 = a + b;
            s2 = b - ( s1 - a );
        end

        function [ s1, s2 ] = TwoSum( a, b )
            if QuadDouble.SingletonExpansionNotSupported
                [ a, b ] = QuadDouble.ExpandSingleton( a, b );
            end
            s1 = a + b;
            bb = s1 - a;
            s2 = ( a - ( s1 - bb ) ) + ( b - bb );
        end

        function [ a, b, c ] = ThreeSum( a, b, c )
            [ t1, t2 ] = QuadDouble.TwoSum( a, b );
            [ a, t3 ]  = QuadDouble.TwoSum( c, t1 );
            [ b, c ]   = QuadDouble.TwoSum( t2, t3 );
        end

        function [ a, b, c ] = ThreeSum2( a, b, c )
            [ t1, t2 ] = QuadDouble.TwoSum( a, b );
            [ a, t3 ]  = QuadDouble.TwoSum( c, t1 );
            b = t2 + t3;
        end

        function [ a1, a2 ] = Split( a )
            t = a .* 134217729; % 2^27 + 1
            a1 = t - ( t - a );
            a2 = a - a1;
        end

        function [ p1, p2 ] = TwoProd( a, b )
            if QuadDouble.SingletonExpansionNotSupported
                [ a, b ] = QuadDouble.ExpandSingleton( a, b );
            end
            p1 = a .* b;
            [ a1, a2 ] = QuadDouble.Split( a );
            [ b1, b2 ] = QuadDouble.Split( b );
            p2 = ( ( a1 .* b1 - p1 ) + a1 .* b2 + a2 .* b1 ) + a2 .* b2;
        end

        function [ p1, p2 ] = TwoSqr( a )
            p1 = a .* a;
            [ a1, a2 ] = QuadDouble.Split( a );
            p2 = ( ( a1 .* a1 - p1 ) + 2 .* a1 .* a2 ) + a2 .* a2;
        end

        function [ c0, c1, c2, c3 ] = Renorm4( c0, c1, c2, c3 )
            [ c3, c4 ] = QuadDouble.QuickTwoSum( c2, c3 );
            [ c2, c3 ] = QuadDouble.QuickTwoSum( c1, c3 );
            [ c1, c2 ] = QuadDouble.QuickTwoSum( c0, c2 );
            c0 = c1;

            s0 = c0;
            s1 = c2;

            s2 = zeros( size( s0 ) );
            s3 = zeros( size( s0 ) );

            Select = s1 ~= 0;
            if any( Select, 'all' )
                [ s1( Select ), s2( Select ) ] = QuadDouble.QuickTwoSum( s1( Select ), c3( Select ) );
                Select2 = Select & ( s2 ~= 0 );
                if any( Select2, 'all' )
                    [ s2( Select2 ), s3( Select2 ) ] = QuadDouble.QuickTwoSum( s2( Select2 ), c4( Select2 ) );
                end
                Select3 = Select & ~( s2 ~= 0 );
                if any( Select3, 'all' )
                    [ s1( Select3 ), s2( Select3 ) ] = QuadDouble.QuickTwoSum( s1( Select3 ), c4( Select3 ) );
                end
            end

            Select4 = ~Select;
            if any( Select4, 'all' )
                [ s0( Select4 ), s1( Select4 ) ] = QuadDouble.QuickTwoSum( s0( Select4 ), c3( Select4 ) );
                Select5 = Select4 & ( s1 ~= 0 );
                if any( Select5, 'all' )
                    [ s1( Select5 ), s2( Select5 ) ] = QuadDouble.QuickTwoSum( s1( Select5 ), c4( Select5 ) );
                end
                Select6 = Select4 & ~( s1 ~= 0 );
                if any( Select6, 'all' )
                    [ s0( Select6 ), s1( Select6 ) ] = QuadDouble.QuickTwoSum( s0( Select6 ), c4( Select6 ) );
                end
            end

            c0 = s0; c1 = s1; c2 = s2; c3 = s3;
        end

        function [ c0, c1, c2, c3 ] = Renorm5( c0, c1, c2, c3, c4 )
            [ c4, c5 ] = QuadDouble.QuickTwoSum( c3, c4 );
            [ c3, c4 ] = QuadDouble.QuickTwoSum( c2, c4 );
            [ c2, c3 ] = QuadDouble.QuickTwoSum( c1, c3 );
            [ c1, c2 ] = QuadDouble.QuickTwoSum( c0, c2 );
            c0 = c1;

            s0 = c0;
            s1 = c2;


            s2 = zeros( size( s0 ) );
            s3 = zeros( size( s0 ) );

            Select = s1 ~= 0;
            if any( Select, 'all' )
                [ s1( Select ), s2( Select ) ] = QuadDouble.QuickTwoSum( s1( Select ), c3( Select ) );
                Select2 = Select & ( s2 ~= 0 );
                if any( Select2, 'all' )
                    [ s2( Select2 ), s3( Select2 ) ] = QuadDouble.QuickTwoSum( s2( Select2 ), c4( Select2 ) );
                    Select2A = Select2 & ( s3 ~= 0 );
                    if any( Select2A, 'all' )
                        s3( Select2A ) = s3( Select2A ) + c5( Select2A );
                    end
                    Select2B = Select2 & ~( s3 ~= 0 );
                    if any( Select2B, 'all' )
                        s2( Select2B ) = s2( Select2B ) + c5( Select2B );
                    end
                end
                Select3 = Select & ~( s2 ~= 0 );
                if any( Select3, 'all' )
                    [ s1( Select3 ), s2( Select3 ) ] = QuadDouble.QuickTwoSum( s1( Select3 ), c4( Select3 ) );
                    Select3A = Select3 & ( s2 ~= 0 );
                    if any( Select3A, 'all' )
                        [ s2( Select3A ), s3( Select3A ) ] = QuadDouble.QuickTwoSum( s2( Select3A ), c5( Select3A ) );
                    end
                    Select3B = Select3 & ~( s2 ~= 0 );
                    if any( Select3B, 'all' )
                        [ s1( Select3B ), s2( Select3B ) ] = QuadDouble.QuickTwoSum( s1( Select3B ), c5( Select3B ) );
                    end
                end
            end

            Select4 = ~Select;
            if any( Select4, 'all' )
                [ s0( Select4 ), s1( Select4 ) ] = QuadDouble.QuickTwoSum( s0( Select4 ), c3( Select4 ) );
                Select5 = Select4 & ( s1 ~= 0 );
                if any( Select5, 'all' )
                    [ s1( Select5 ), s2( Select5 ) ] = QuadDouble.QuickTwoSum( s1( Select5 ), c4( Select5 ) );
                    Select5A = Select5 & ( s2 ~= 0 );
                    if any( Select5A, 'all' )
                        [ s2( Select5A ), s3( Select5A ) ] = QuadDouble.QuickTwoSum( s2( Select5A ), c5( Select5A ) );
                    end
                    Select5B = Select5 & ~( s2 ~= 0 );
                    if any( Select5B, 'all' )
                        [ s1( Select5B ), s2( Select5B ) ] = QuadDouble.QuickTwoSum( s1( Select5B ), c5( Select5B ) );
                    end
                end
                Select6 = Select4 & ~( s1 ~= 0 );
                if any( Select6, 'all' )
                    [ s0( Select6 ), s1( Select6 ) ] = QuadDouble.QuickTwoSum( s0( Select6 ), c4( Select6 ) );
                    Select6A = Select6 & ( s1 ~= 0 );
                    if any( Select6A, 'all' )
                        [ s1( Select6A ), s2( Select6A ) ] = QuadDouble.QuickTwoSum( s1( Select6A ), c5( Select6A ) );
                    end
                    Select6B = Select6 & ~( s1 ~= 0 );
                    if any( Select6B, 'all' )
                        [ s0( Select6B ), s1( Select6B ) ] = QuadDouble.QuickTwoSum( s0( Select6B ), c5( Select6B ) );
                    end
                end
            end

            c0 = s0; c1 = s1; c2 = s2; c3 = s3;
        end

        function [ s0, s1, s2, s3 ] = QDPlusQD( a0, a1, a2, a3, b0, b1, b2, b3 )
            if QuadDouble.SingletonExpansionNotSupported
                [ a0, b0 ] = QuadDouble.ExpandSingleton( a0, b0 );
                [ a1, b1 ] = QuadDouble.ExpandSingleton( a1, b1 );
                [ a2, b2 ] = QuadDouble.ExpandSingleton( a2, b2 );
                [ a3, b3 ] = QuadDouble.ExpandSingleton( a3, b3 );
            end
            ss0 = a0 + b0;
            ss1 = a1 + b1;
            ss2 = a2 + b2;
            ss3 = a3 + b3;
            v0 = ss0 - a0;
            v1 = ss1 - a1;
            v2 = ss2 - a2;
            v3 = ss3 - a3;
            u0 = ss0 - v0;
            u1 = ss1 - v1;
            u2 = ss2 - v2;
            u3 = ss3 - v3;
            w0 = a0 - u0;
            w1 = a1 - u1;
            w2 = a2 - u2;
            w3 = a3 - u3;
            uu0 = b0 - v0;
            uu1 = b1 - v1;
            uu2 = b2 - v2;
            uu3 = b3 - v3;
            t0 = w0 + uu0;
            t1 = w1 + uu1;
            t2 = w2 + uu2;
            t3 = w3 + uu3;
            [ ss1, t0_2 ] = QuadDouble.TwoSum( ss1, t0 );
            [ ss2, t0_3, t1_2 ] = QuadDouble.ThreeSum( ss2, t0_2, t1 );
            [ ss3, t0_4, ~ ] = QuadDouble.ThreeSum2( ss3, t0_3, t2 );
            t0_5 = t0_4 + t1_2 + t3;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( ss0, ss1, ss2, ss3, t0_5 );
        end

        function [ s0, s1, s2, s3 ] = QDPlusDouble( a0, a1, a2, a3, b )
            if QuadDouble.SingletonExpansionNotSupported
                [ a0, b ] = QuadDouble.ExpandSingleton( a0, b );
                [ a1, b ] = QuadDouble.ExpandSingleton( a1, b );
                [ a2, b ] = QuadDouble.ExpandSingleton( a2, b );
                [ a3, b ] = QuadDouble.ExpandSingleton( a3, b );
            end
            [ c0, e ] = QuadDouble.TwoSum( a0, b );
            [ c1, e ] = QuadDouble.TwoSum( a1, e );
            [ c2, e ] = QuadDouble.TwoSum( a2, e );
            [ c3, ~ ] = QuadDouble.TwoSum( a3, e );
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( c0, c1, c2, c3, e );
        end

        function [ s0, s1, s2, s3 ] = QDTimesQD( a0, a1, a2, a3, b0, b1, b2, b3 )
            if QuadDouble.SingletonExpansionNotSupported
                [ a0, b0 ] = QuadDouble.ExpandSingleton( a0, b0 );
                [ a1, b1 ] = QuadDouble.ExpandSingleton( a1, b1 );
                [ a2, b2 ] = QuadDouble.ExpandSingleton( a2, b2 );
                [ a3, b3 ] = QuadDouble.ExpandSingleton( a3, b3 );
            end
            [ p0, q0 ] = QuadDouble.TwoProd( a0, b0 );
            [ p1, q1 ] = QuadDouble.TwoProd( a0, b1 );
            [ p2, q2 ] = QuadDouble.TwoProd( a1, b0 );
            [ p3, q3 ] = QuadDouble.TwoProd( a0, b2 );
            [ p4, q4 ] = QuadDouble.TwoProd( a1, b1 );
            [ p5, q5 ] = QuadDouble.TwoProd( a2, b0 );
            [ p1, p2, q0 ] = QuadDouble.ThreeSum( p1, p2, q0 );
            [ p2, q1, q2 ] = QuadDouble.ThreeSum( p2, q1, q2 );
            [ p3, p4, p5 ] = QuadDouble.ThreeSum( p3, p4, p5 );
            [ ss0, t0 ] = QuadDouble.TwoSum( p2, p3 );
            [ ss1, t1 ] = QuadDouble.TwoSum( q1, p4 );
            ss2 = q2 + p5;
            [ ss1, t0_2 ] = QuadDouble.TwoSum( ss1, t0 );
            ss2 = ss2 + ( t0_2 + t1 );
            ss1 = ss1 + a0.*b3 + a1.*b2 + a2.*b1 + a3.*b0 + q0 + q3 + q4 + q5;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( p0, p1, ss0, ss1, ss2 );
        end

        function [ s0, s1, s2, s3 ] = QDTimesDouble( a0, a1, a2, a3, b )
            if QuadDouble.SingletonExpansionNotSupported
                [ a0, b ] = QuadDouble.ExpandSingleton( a0, b );
                [ a1, b ] = QuadDouble.ExpandSingleton( a1, b );
                [ a2, b ] = QuadDouble.ExpandSingleton( a2, b );
                [ a3, b ] = QuadDouble.ExpandSingleton( a3, b );
            end
            [ p0, q0 ] = QuadDouble.TwoProd( a0, b );
            [ p1, q1 ] = QuadDouble.TwoProd( a1, b );
            [ p2, q2 ] = QuadDouble.TwoProd( a2, b );
            p3 = a3 .* b;
            ss0 = p0;
            [ ss1, s2_1 ] = QuadDouble.TwoSum( q0, p1 );
            [ s2_1, q1, p2 ] = QuadDouble.ThreeSum( s2_1, q1, p2 );
            [ q1, q2, ~ ] = QuadDouble.ThreeSum2( q1, q2, p3 );
            ss3 = q1;
            s4 = q2 + p2;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( ss0, ss1, s2_1, ss3, s4 );
        end

        function [ s0, s1, s2, s3 ] = QDDivQD( a0, a1, a2, a3, b0, b1, b2, b3 )
            q0 = a0 ./ b0;
            [ r0, r1, r2, r3 ] = QuadDouble.QDTimesDouble( b0, b1, b2, b3, q0 );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( a0, a1, a2, a3, -r0, -r1, -r2, -r3 );
            q1 = r0 ./ b0;
            [ t0, t1, t2, t3 ] = QuadDouble.QDTimesDouble( b0, b1, b2, b3, q1 );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
            q2 = r0 ./ b0;
            [ t0, t1, t2, t3 ] = QuadDouble.QDTimesDouble( b0, b1, b2, b3, q2 );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
            q3 = r0 ./ b0;
            [ t0, t1, t2, t3 ] = QuadDouble.QDTimesDouble( b0, b1, b2, b3, q3 );
            [ r0, ~, ~, ~ ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -t0, -t1, -t2, -t3 );
            q4 = r0 ./ b0;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( q0, q1, q2, q3, q4 );
        end

        function [ s0, s1, s2, s3 ] = QDDivDouble( a0, a1, a2, a3, b )
            q0 = a0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q0, b );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( a0, a1, a2, a3, -p0, -p1, 0, 0 );
            q1 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q1, b );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -p0, -p1, 0, 0 );
            q2 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q2, b );
            [ r0, r1, r2, r3 ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -p0, -p1, 0, 0 );
            q3 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q3, b );
            [ r0, ~, ~, ~ ] = QuadDouble.QDPlusQD( r0, r1, r2, r3, -p0, -p1, 0, 0 );
            q4 = r0 ./ b;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( q0, q1, q2, q3, q4 );
        end

        function [ s0, s1, s2, s3 ] = DoubleDivDouble( a, b )
            q0 = a ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q0, b );
            [ r0, ~ ] = QuadDouble.TwoSum( a, -p0 );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p1 );
            q1 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q1, b );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p0 );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p1 );
            q2 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q2, b );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p0 );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p1 );
            q3 = r0 ./ b;
            [ p0, p1 ] = QuadDouble.TwoProd( q3, b );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p0 );
            [ r0, ~ ] = QuadDouble.TwoSum( r0, -p1 );
            q4 = r0 ./ b;
            [ s0, s1, s2, s3 ] = QuadDouble.Renorm5( q0, q1, q2, q3, q4 );
        end

        function v = Index( v, s )
            if isa( v, 'QuadDouble' )
                v = QuadDouble.Make( v.v1.v1( s ), v.v1.v2( s ), v.v2.v1( s ), v.v2.v2( s ) );
            else
                v = v( s );
            end
        end

        function v = Assign( v, s, x )
            if isa( v, 'QuadDouble' ) || isa( x, 'QuadDouble' )
                v = QuadDouble( v );
                x = QuadDouble( x );
                v.v1.v1( s ) = x.v1.v1;
                v.v1.v2( s ) = x.v1.v2;
                v.v2.v1( s ) = x.v2.v1;
                v.v2.v2( s ) = x.v2.v2;
            else
                v( s ) = x;
            end
        end

        function Supported = TestSingletonExpansion
            a = ones( 2, 1 );
            b = [];
            c = [];
            try
                b = a + a.';
                c = a .* a.';
            catch
            end
            Supported = all( [ size( b ), size( c ) ] == 2 );
        end

    end
end
