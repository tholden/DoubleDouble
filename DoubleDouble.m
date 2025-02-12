% The below code is derived from the QD C++ library.

% QD is Copyright (c) 2003-2009, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from U.S. Dept. of Energy) All rights reserved.

% QD is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% (1) Redistributions of source code must retain the copyright notice, this list of conditions and the following disclaimer.
% (2) Redistributions in binary form must reproduce the copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 3. You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.

% The implementation of the LUP decomposition and the backslack operator here is derived from Cleve Moler's code from "Numerical Computing with MATLAB".

% "Numerical Computing with MATLAB" is Copyright (c) 2004, Cleve Moler and Copyright (c) 2016, The MathWorks, Inc.

% "Numerical Computing with MATLAB" is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% (3) In all cases, the software is, and all modifications and derivatives of the software shall be, licensed to you solely for use in conjunction with MathWorks products and service offerings.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

% The implementation of the LDL decomposition here is derived from Brian Borchers's ldlt library.

% Brian Borchers's ldlt library is Copyright (c) 2009, Brian Borchers.

% Brian Borchers's ldlt library is distributed under the following license:

% 1. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
% (3) Neither the name of the New Mexico Inst of Mining & Tech nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

classdef DoubleDouble
    properties ( SetAccess = private, GetAccess = private )
        v1
        v2
    end
    
    properties ( Constant, GetAccess = public )
        zero = DoubleDouble.Make( 0, 0 );
        one = DoubleDouble.Make( 1, 0 );
        eps = DoubleDouble.Make( 4.93038065763132e-32, 0 );
        pi = DoubleDouble.Make( 3.141592653589793116e+00, 1.224646799147353207e-16 );
    end
    
    properties ( Constant, GetAccess = private )
        SingletonExpansionNotSupported = ~DoubleDouble.TestSingletonExpansion( );
        
        InverseFactorial = [
            1.66666666666666657e-01,  9.25185853854297066e-18;
            4.16666666666666644e-02,  2.31296463463574266e-18;
            8.33333333333333322e-03,  1.15648231731787138e-19;
            1.38888888888888894e-03, -5.30054395437357706e-20;
            1.98412698412698413e-04,  1.72095582934207053e-22;
            2.48015873015873016e-05,  2.15119478667758816e-23;
            2.75573192239858925e-06, -1.85839327404647208e-22;
            2.75573192239858883e-07,  2.37677146222502973e-23;
            2.50521083854417202e-08, -1.44881407093591197e-24;
            2.08767569878681002e-09, -1.20734505911325997e-25;
            1.60590438368216133e-10,  1.25852945887520981e-26;
            1.14707455977297245e-11,  2.06555127528307454e-28;
            7.64716373181981641e-13,  7.03872877733453001e-30;
            4.77947733238738525e-14,  4.39920548583408126e-31;
            2.81145725434552060e-15,  1.65088427308614326e-31;
        ];
    
        NInverseFactorial = 15;
        
        piT2   = DoubleDouble.Make(  6.283185307179586232e+00,  2.449293598294706414e-16 );
        piD2   = DoubleDouble.Make(  1.570796326794896558e+00,  6.123233995736766036e-17 );
        piD16  = DoubleDouble.Make(  1.963495408493620697e-01,  7.654042494670957545e-18 );
        
        log_2  = DoubleDouble.Make(  6.931471805599452862e-01,  2.319046813846299558e-17 );
        log_10 = DoubleDouble.Make(  2.302585092994045901e+00, -2.170756223382249351e-16 );
        
        SinTable = [
            0, 0;
            1.950903220161282758e-01, -7.991079068461731263e-18;
            3.826834323650897818e-01, -1.005077269646158761e-17;
            5.555702330196021776e-01,  4.709410940561676821e-17;
            7.071067811865475727e-01, -4.833646656726456726e-17;
        ];
    
        CosTable = [
            0, 0;
            9.807852804032304306e-01,  1.854693999782500573e-17;
            9.238795325112867385e-01,  1.764504708433667706e-17;
            8.314696123025452357e-01,  1.407385698472802389e-18;
            7.071067811865475727e-01, -4.833646656726456726e-17;
        ];
    end
    
    methods
        function v = DoubleDouble( in, varargin )
            if nargin == 0
                v.v1 = [];
                v.v2 = [];
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = DoubleDouble.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'DoubleDouble' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            else
                v.v1 = double( in );
                v.v2 = zeros( size( in ) );
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
        
        function [ x1, x2 ] = ToSumOfDoubles( v )
            x1 = v.v1;
            x2 = v.v2;
        end
        
        function v = double( v )
            v = v.v1;
        end
        
        function v = isreal( v )
            v = isreal( v.v1 ) && isreal( v.v2 );
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
                v.v1( Select ) = 0;
                v.v2( Select ) = 0;
                v.v1( ~Select ) = DoubleDouble.pi.v1;
                v.v2( ~Select ) = DoubleDouble.pi.v2;
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
            n = numArgumentsFromSubscript( v.v1, s, IndexingContext );
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
            v = DoubleDouble.Make( repmat( v.v1, varargin{:} ), repmat( v.v2, varargin{:} ) );
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
            v = isempty( v.v1 );
        end
        
        function v = diag( v, k )
            if nargin < 2
                v = DoubleDouble.Make( diag( v.v1 ), diag( v.v2 ) );
            else
                v = DoubleDouble.Make( diag( v.v1, k ), diag( v.v2, k ) );
            end
        end
        
        function v = tril( v, k )
            if nargin < 2
                v = DoubleDouble.Make( tril( v.v1 ), tril( v.v2 ) );
            else
                v = DoubleDouble.Make( tril( v.v1, k ), tril( v.v2, k ) );
            end
        end
        
        function v = triu( v, k )
            if nargin < 2
                v = DoubleDouble.Make( triu( v.v1 ), triu( v.v2 ) );
            else
                v = DoubleDouble.Make( triu( v.v1, k ), triu( v.v2, k ) );
            end
        end
        
        function v = plus( a, b )
            v = DoubleDouble.Plus( a, b );
        end
        
        function v = minus( a, b )
            v = DoubleDouble.Minus( a, b );
        end
        
        function v = uminus( v )
            v.v1 = -v.v1;
            v.v2 = -v.v2;
        end
        
        function v = uplus( v )
        end
        
        function v = times( a, b )
            v = DoubleDouble.Times( a, b );
        end
                
        function v = mtimes( a, b )
            v = DoubleDouble.MTimes( a, b );
        end
        
        function v = rdivide( a, b )
            v = DoubleDouble.RDivide( a, b );
        end
                
        function v = ldivide( a, b )
            v = DoubleDouble.LDivide( a, b );
        end
        
        function v = mldivide( a, v )
            v = DoubleDouble.MLDivide( a, v );
        end
        
        function v = mrdivide( v, a )
            v = DoubleDouble.MRDivide( v, a );
        end
        
        function v = power( a, b )
            if ~isa( a, 'DoubleDouble' )
                a = DoubleDouble( a );
            end
            v = exp( b .* log( a ) );
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
                    v = DoubleDouble;
                end
            end
        end
        
        function v = lt( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( a2 < b2 ) );
        end
        
        function v = gt( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( a2 > b2 ) );
        end
        
        function v = le( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 < b1 ) | ( ( a1 == b1 ) & ( a2 <= b2 ) );
        end
        
        function v = ge( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 > b1 ) | ( ( a1 == b1 ) & ( a2 >= b2 ) );
        end
        
        function v = ne( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
                b1 = b.v1;
                b2 = b.v2;
            else
                b1 = b;
                b2 = 0;
            end
            v = ( a1 ~= b1 ) | ( a2 ~= b2 );
        end
        
        function v = eq( a, b )
            if isa( a, 'DoubleDouble' )
                a1 = a.v1;
                a2 = a.v2;
            else
                a1 = a;
                a2 = 0;
            end
            if isa( b, 'DoubleDouble' )
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
            if ~isa( a, 'DoubleDouble' )
                a = DoubleDouble( a );
            end
            if ~isa( b, 'DoubleDouble' )
                b = DoubleDouble( b );
            end
            if ~isa( d, 'DoubleDouble' )
                d = DoubleDouble( d );
            end
            c = double( floor( ( b - a ) ./ d ) );
            v = a + ( 0:c ) .* d;
        end
        
        function v = ctranspose( v )
            v.v1 = v.v1';
            v.v2 = v.v2';
        end
        
        function v = transpose( v )
            v.v1 = v.v1.';
            v.v2 = v.v2.';
        end
        
        function v = horzcat( a, b, varargin )
            if nargin > 2
                v = horzcat( horzcat( a, b ), varargin{:} );
            else
                if ~isa( a, 'DoubleDouble' )
                    a = DoubleDouble( a );
                end
                if ~isa( b, 'DoubleDouble' )
                    b = DoubleDouble( b );
                end
                x1 = horzcat( [ a.v1, b.v1 ] );
                x2 = horzcat( [ a.v2, b.v2 ] );
                v = DoubleDouble.Make( x1, x2 );
            end
        end
        
        function v = vertcat( a, b, varargin )
            if nargin > 2
                v = vertcat( vertcat( a, b ), varargin{:} );
            else
                if ~isa( a, 'DoubleDouble' )
                    a = DoubleDouble( a );
                end
                if ~isa( b, 'DoubleDouble' )
                    b = DoubleDouble( b );
                end
                x1 = vertcat( [ a.v1; b.v1 ] );
                x2 = vertcat( [ a.v2; b.v2 ] );
                v = DoubleDouble.Make( x1, x2 );
            end
        end
        
        function v = subsref( v, s )
            v.v1 = subsref( v.v1, s );
            v.v2 = subsref( v.v2, s );
        end
        
        function v = subsasgn( v, s, b )
            if ~isa( v, 'DoubleDouble' )
                v = DoubleDouble( v );
            end
            if ~isa( b, 'DoubleDouble' )
                b = DoubleDouble( b );
            end
            v.v1 = subsasgn( v.v1, s, b.v1 );
            v.v2 = subsasgn( v.v2, s, b.v2 );
        end
        
        function v = subsindex( v )
            v = v.v1;
        end
        
        function [ v, Indices ] = sort( v, varargin )
            DimIndex = find( cellfun( @isnumeric, varargin ), 1 );
            if isempty( DimIndex )
                dim = [];
            else
                dim = varargin{ DimIndex };
                varargin = varargin( [ 1 : ( DimIndex - 1 ), ( DimIndex + 1 ) : end ] );
            end
            CMIndex = find( strcmpi( varargin, 'ComparisonMethod' ), 1 );
            if isempty( CMIndex )
                cm = [];
            else
                cm = varargin{ CMIndex + 1 };
                varargin = varargin( [ 1 : ( CMIndex - 1 ), ( CMIndex + 2 ) : end ] );
            end
            [ v, Indices ] = DoubleDouble.Sort( v, dim, cm, varargin{:} );
        end
        
        function v = sum( v, dim )
            if nargin < 2
                dim = [];
            end
            v = DoubleDouble.Sum( v, dim );
        end
        
        function v = prod( v, dim )
            if nargin < 2
                dim = [];
            end
            v = DoubleDouble.Prod( v, dim );
        end
        
        function [ v, i ] = max( a, b, dim )
            if nargin < 3
                dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if nargout < 2
                v = DoubleDouble.Max( a, b, dim );
            else
                [ v, i ] = DoubleDouble.Max( a, b, dim );
            end
        end
        
        function [ v, i ] = min( a, b, dim )
            if nargin < 3
                dim = [];
                if nargin < 2
                    b = [];
                end
            end
            if nargout < 2
                v = DoubleDouble.Min( a, b, dim );
            else
                [ v, i ] = DoubleDouble.Min( a, b, dim );
            end
        end
        
        function v = cumsum( v, dim )
            if nargin < 2
                dim = [];
            end
            v = DoubleDouble.CumSum( v, dim );
        end
        
        function v = diff( v, dim )
            if nargin < 2
                dim = [];
            end
            v = DoubleDouble.Diff( v, dim );
        end
        
        function v = cumprod( v, dim )
            if nargin < 2
                dim = [];
            end
            v = DoubleDouble.CumProd( v, dim );
        end
        
        function v = cummax( v, dim )
            if nargin < 3
                dim = [];
            end
            v = DoubleDouble.CumMax( v, dim );
        end
        
        function v = cummin( v, dim )
            if nargin < 3
                dim = [];
            end
            v = DoubleDouble.CumMin( v, dim );
        end
        
        function v = dot( a, b, dim )
            if nargin < 3
                dim = [];
            end
            v = DoubleDouble.Dot( a, b, dim );
        end
        
        function v = norm( v, p )
            if nargin < 2
                v = DoubleDouble.Norm( v );
            else
                v = DoubleDouble.Norm( v, p );
            end
        end
        
        function v = abs( v )
            if isreal( v )
                Select = v.v1 < 0;
                v.v1( Select ) = -v.v1( Select );
                v.v2( Select ) = -v.v2( Select ); 
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
                [ v.v1, v.v2 ] = DoubleDouble.DDDividedByDD( v.v1, v.v2, abs_v.v1, abs_v.v2, true );
            end
        end
        
        function v = floor( v )
            x1 = floor( v.v1 );
            x2 = zeros( size( x1 ) );
            Select = x1 == v.v1;
            x2( Select ) = floor( v.v2( Select ) );
            [ x1, x2 ] = DoubleDouble.Normalize( x1, x2 );
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = ceil( v )
            x1 = ceil( v.v1 );
            x2 = zeros( size( x1 ) );
            Select = x1 == v.v1;
            x2( Select ) = ceil( v.v2( Select ) );
            [ x1, x2 ] = DoubleDouble.Normalize( x1, x2 );
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = fix( v )
            x1 = fix( v.v1 );
            x2 = zeros( size( x1 ) );
            Select = x1 == v.v1;
            x2( Select ) = fix( v.v2( Select ) );
            [ x1, x2 ] = DoubleDouble.Normalize( x1, x2 );
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = round( v )
            x1 = round( v.v1 );
            x2 = zeros( size( x1 ) );
            Select = x1 == v.v1;
            x2( Select ) = round( v.v2( Select ) );
            Select = ( ~Select ) & ( abs( x1 - v.v1 ) == 0.5 ) & ( v.v2 < 0 );
            x2( Select ) = x2( Select ) - 1;
            [ x1, x2 ] = DoubleDouble.Normalize( x1, x2 );
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = realsqrt( v )
            Select = v < 0;
            v.v1( Select ) = NaN;
            v.v2( Select ) = NaN;
            Select = v > 0;
            x = 1 ./ sqrt( v.v1( Select ) );
            vx = v.v1( Select ) .* x;
            t = DoubleDouble.Make( v.v1( Select ), v.v2( Select ) ) - DoubleDouble.Times( vx, vx );
            t = DoubleDouble.Plus( vx, t.v1 .* ( x * 0.5 ) );
            v.v1( Select ) = t.v1;
            v.v2( Select ) = t.v2;
        end
        
        function v = sqrt( v )
            Select = v ~= 0;
            x = 1 ./ sqrt( v.v1( Select ) );
            vx = v.v1( Select ) .* x;
            t = DoubleDouble.Make( v.v1( Select ), v.v2( Select ) ) - DoubleDouble.Times( vx, vx );
            t = DoubleDouble.Plus( vx, t.v1 .* ( x * 0.5 ) );
            v.v1( Select ) = t.v1;
            v.v2( Select ) = t.v2;
        end
        
        function v = sqrtm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( sqrt( d ) ) / v;
        end
        
        function v = exp( v )
            if ~isreal( v )
                [ sin_imag_v, cos_imag_v ] = sincos( imag( v ) );
                v = exp( real( v ) ) .* ( cos_imag_v + 1i .* sin_imag_v );
                return
            end
            
            % Strategy:  We first reduce the size of x by noting that
            % exp(kr + m * log(2)) = 2^m * exp(r)^k
            % where m and k are integers.  By choosing m appropriately
            % we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is 
            % evaluated using the familiar Taylor series.  Reducing the 
            % argument substantially speeds up the convergence.
            k = 512.0;
            inv_k = 1.0 / k;
            Threshhold = inv_k .* DoubleDouble.eps.v1;

            m = floor( v.v1 ./ DoubleDouble.log_2.v1 + 0.5 );
            r = TimesPowerOf2( v - DoubleDouble.log_2 .* m, inv_k );

            p = r .* r;
            s = r + TimesPowerOf2( p, 0.5 );
            p = p .* r;
            t = p .* DoubleDouble.Make( DoubleDouble.InverseFactorial( 1, 1 ), DoubleDouble.InverseFactorial( 1, 2 ) );
            for i = 2 : DoubleDouble.NInverseFactorial
                s = s + t;
                p = p .* r;
                t = p .* DoubleDouble.Make( DoubleDouble.InverseFactorial( i, 1 ), DoubleDouble.InverseFactorial( i, 2 ) );
                if all( abs( t.v1(:) ) <= Threshhold )
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
            s = s + 1.0;
            
            v = DoubleDouble.Make( pow2( s.v1, m ), pow2( s.v2, m ) );
        end
        
        function v = expm( v )
            [ v, d ] = eig( v );
            d = diag( d );
            assert( length( unique( d.v1 ) ) == length( d.v1 ) );
            v = v * diag( exp( d ) ) / v;
        end
        
        function x = log( v )
            x = DoubleDouble.Make( log( v.v1 ), zeros( size( v.v1 ) ) );
            x = x + v .* exp( -x ) - 1.0;
            x = x + v .* exp( -x ) - 1.0; % slightly paranoid, but does correct e.g. log(exp(DoubleDouble(-40)))
        end
        
        function v = log2( v )
            v = log( v ) ./ DoubleDouble.log_2;
        end
        
        function v = log10( v )
            v = log( v ) ./ DoubleDouble.log_10;
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
                exp_Piv = exp( TimesPowerOf2( v, 1i ) );
                exp_Niv = 1 ./ exp_Piv;
                sin_v = TimesPowerOf2( exp_Piv - exp_Niv, -0.5i );
                cos_v = TimesPowerOf2( exp_Piv + exp_Niv, +0.5  );
                return
            end
            % Strategy.  To compute sin(x), cos(x), we choose integers a, b so that
            % x = s + a * (pi/2) + b * (pi/16)
            % and |s| <= pi/32.  Using the fact that 
            % sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))
            % we can compute sin(x) from sin(s), cos(s).  This greatly increases the convergence of the sine Taylor series.
            
            z = round( v ./ DoubleDouble.piT2 );
            r = v - DoubleDouble.piT2 .* z;
            
            q = floor( r.v1 ./ DoubleDouble.piD2.v1 + 0.5 );
            t = r - DoubleDouble.piD2 .* q;
            j = q;
            abs_j = abs( j );
            
            q = floor( t.v1 ./ DoubleDouble.piD16.v1 + 0.5 );
            t = t - DoubleDouble.piD16 .* q;
            k = q;
            abs_k = abs( k );
            
            test = ( j >= -2 ) & ( j <= 2 );
            assert( all( test(:) ) );
            test = abs_k <= 4;
            assert( all( test(:) ) );
            
            [ sin_t, cos_t ] = SinCosTaylor( t );
            
            sin_v = sin_t;
            cos_v = cos_t;
            
            a = DoubleDouble.Make( DoubleDouble.CosTable( abs_k + 1, 1 ), DoubleDouble.CosTable( abs_k + 1, 2 ) );
            b = DoubleDouble.Make( DoubleDouble.SinTable( abs_k + 1, 1 ), DoubleDouble.SinTable( abs_k + 1, 2 ) );
            
            a = reshape( a, size( v ) );
            b = reshape( b, size( v ) );

            a_sin_t = a .* sin_t;
            b_sin_t = b .* sin_t;
            a_cos_t = a .* cos_t;
            b_cos_t = b .* cos_t;
            
            Select = k > 0;
            
            [ sin_v.v1( Select ), sin_v.v2( Select ) ] = DoubleDouble.DDPlusDD( +a_sin_t.v1( Select ), +a_sin_t.v2( Select ), +b_cos_t.v1( Select ), +b_cos_t.v2( Select ) );
            [ cos_v.v1( Select ), cos_v.v2( Select ) ] = DoubleDouble.DDPlusDD( -b_sin_t.v1( Select ), -b_sin_t.v2( Select ), +a_cos_t.v1( Select ), +a_cos_t.v2( Select ) );
            
            Select = k < 0;
            
            [ sin_v.v1( Select ), sin_v.v2( Select ) ] = DoubleDouble.DDPlusDD( +a_sin_t.v1( Select ), +a_sin_t.v2( Select ), -b_cos_t.v1( Select ), -b_cos_t.v2( Select ) );
            [ cos_v.v1( Select ), cos_v.v2( Select ) ] = DoubleDouble.DDPlusDD( +b_sin_t.v1( Select ), +b_sin_t.v2( Select ), +a_cos_t.v1( Select ), +a_cos_t.v2( Select ) );

            Select = j == 1;
            
            [ sin_v.v1( Select ), sin_v.v2( Select ), cos_v.v1( Select ), cos_v.v2( Select ) ] = deal( +cos_v.v1( Select ), +cos_v.v2( Select ), -sin_v.v1( Select ), -sin_v.v2( Select ) );
            
            Select = j == -1;

            [ sin_v.v1( Select ), sin_v.v2( Select ), cos_v.v1( Select ), cos_v.v2( Select ) ] = deal( -cos_v.v1( Select ), -cos_v.v2( Select ), +sin_v.v1( Select ), +sin_v.v2( Select ) );
            
            Select = abs_j == 2;
            
            [ sin_v.v1( Select ), sin_v.v2( Select ), cos_v.v1( Select ), cos_v.v2( Select ) ] = deal( -sin_v.v1( Select ), -sin_v.v2( Select ), -cos_v.v1( Select ), -cos_v.v2( Select ) );
        end
        
        function v = sin( v )
            [ v, ~ ] = sincos( v );
        end
        
        function v = asin( v )
            assert( all( abs( v.v1(:) ) <= 1 ) );
            v = atan2( v, sqrt( 1 - v.*v ) );
        end
        
        function v = cos( v )
            [ ~, v ] = sincos( v );
        end
        
        function v = acos( v )
            assert( all( abs( v.v1(:) ) <= 1 ) );
            v = atan2( sqrt( 1 - v.*v ), v );
        end
        
        function v = tan( v )
            [ sin_v, cos_v ] = sincos( v );
            v = sin_v ./ cos_v;
        end
        
        function v = atan( v )
            v = atan2( v, DoubleDouble.Make( 1, 0 ) );
        end
        
        function v = atan2( y, x )
            r = sqrt( x.*x + y.*y );
            xx = x ./ r;
            yy = y ./ r;
            Select = abs( xx.v1 ) > abs( yy.v1 );
            v = DoubleDouble( atan2( y.v1, x.v1 ) );
            [ sin_z, cos_z ] = sincos( v );
            t = yy;
            [ t.v1( Select ), t.v2( Select ) ] = DoubleDouble.DDPlusDD( t.v1( Select ), t.v2( Select ),      -sin_z.v1( Select ), -sin_z.v2( Select ) );
            [ t.v1( Select ), t.v2( Select ) ] = DoubleDouble.DDDividedByDD( t.v1( Select ), t.v2( Select ), +cos_z.v1( Select ), +cos_z.v2( Select ) );
            Select = ~Select;
            [ t.v1( Select ), t.v2( Select ) ] = DoubleDouble.DDPlusDD( xx.v1( Select ), xx.v2( Select ),    -cos_z.v1( Select ), -cos_z.v2( Select ) );
            [ t.v1( Select ), t.v2( Select ) ] = DoubleDouble.DDDividedByDD( t.v1( Select ), t.v2( Select ), -sin_z.v1( Select ), -sin_z.v2( Select ) );
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
                [ ~, midx ] = max( abs( DoubleDouble.Make( v.v1( k:m, k ), v.v2( k:m, k ) ) ) );
                midx = midx + k - 1;

                % Skip elimination if column is zero
                if v.v1( midx, k ) ~= 0 || v.v2( midx, k ) ~= 0

                    % Swap pivot row
                    if midx ~= k
                        v.v1( [ k midx ], : ) = v.v1( [ midx k ], : );
                        v.v2( [ k midx ], : ) = v.v2( [ midx k ], : );
                        p( [ k midx ] ) = p( [ midx k ] );
                    end

                    % Compute multipliers
                    i = k + 1 : m;
                    [ v.v1( i, k ), v.v2( i, k ) ] = DoubleDouble.DDDividedByDD( v.v1( i, k ), v.v2( i, k ), v.v1( k, k ), v.v2( k, k ) );

                    % Update the remainder of the matrix
                    j = k + 1 : n;
                    % A( i, j ) = A( i, j ) - A( i, k ) .* A( k, j );
                    [ t1, t2 ] = DoubleDouble.DDTimesDD( v.v1( i, k ), v.v2( i, k ), v.v1( k, j ), v.v2( k, j ) );
                    [ v.v1( i, j ), v.v2( i, j ) ] = DoubleDouble.DDPlusDD( v.v1( i, j ), v.v2( i, j ), -t1, -t2 );
                end
            end

            if nargout > 1
                % Separate result
                L = tril( v, -1 ) + eye( m, n, 'DoubleDouble' );
                U = triu( v );
                if n > m
                    L.v1 = L.v1( :, 1:m );
                    L.v2 = L.v2( :, 1:m );
                elseif n < m
                    U.v1 = U.v1( 1:n, : );
                    U.v2 = U.v2( 1:n, : );
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
                    v.v1 = v.v1( invp, : );
                    v.v2 = v.v2( invp, : );
                end
            end            
        end
        
        function [ q, v ] = qr( v )
            [ m, n ] = size( v );
            I = eye( m, 'DoubleDouble' );
            QT = I;
            for c = 1 : min( m - 1, n )
                x = DoubleDouble.Make( v.v1( :, c ), v.v2( :, c ) );
                x.v1( 1 : ( c - 1 ) ) = 0;
                x.v2( 1 : ( c - 1 ) ) = 0;
                alpha = norm( x );
                sign_x_c = sign( DoubleDouble.Make( x.v1( c ), x.v2( c ) ) );
                if sign_x_c ~= 0
                    alpha = -alpha .* sign_x_c;
                end
                a = x;
                [ a.v1( c ), a.v2( c ) ] = DoubleDouble.DDPlusDD( a.v1( c ), a.v2( c ), -alpha.v1, -alpha.v2 );
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
                v = DoubleDouble.Make( NaN, NaN );
            end            
        end         
        
        function v = inv( v )
            n = size( v, 1 );
            v = v \ DoubleDouble.Make( eye( n ), zeros( n ) );
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
            L = DoubleDouble.Make( eye( n ), zeros( n ) );
            x1 = zeros( 1, n );
            x2 = x1;
            t1 = x1;
            t2 = x1;
            d = DoubleDouble.Make( x1, x1 );
            x1( 1 ) = v.v1( 1, 1 );
            x2( 1 ) = v.v2( 1, 1 );
            d.v1( 1 ) = x1( 1 );
            d.v2( 1 ) = x2( 1 );
            idxs = 2 : n;
            [ L.v1( idxs, 1 ), L.v2( idxs, 1 ) ] = DoubleDouble.DDDividedByDD( v.v1( idxs, 1 ), v.v2( 2 : n, 1 ), x1( 1 ), x2( 1 ) );
            for j = 2 : n
                idxs = 1 : j - 1;
                [ x1( idxs ), x2( idxs ) ] = DoubleDouble.DDTimesDD( conj( L.v1( j, idxs ) ), conj( L.v2( j, idxs ) ), d.v1( idxs ), d.v2( idxs ) );
                [ t1( idxs ), t2( idxs ) ] = DoubleDouble.DDTimesDD( L.v1( j, idxs ), L.v2( j, idxs ), x1( idxs ), x2( idxs ) );
                t = sum( DoubleDouble.Make( t1( idxs ), t2( idxs ) ) );
                [ x1( j ), x2( j ) ] = DoubleDouble.DDPlusDD( v.v1( j, j ), v.v2( j, j ), -t.v1, -t.v2 );
                d.v1( j ) = x1( j );
                d.v2( j ) = x2( j );
                if j < n
                    jdxs = j + 1 : n;
                    [ s1, s2 ] = DoubleDouble.DDTimesDD( L.v1( jdxs, idxs ), L.v2( jdxs, idxs ), x1( idxs ), x2( idxs ) );
                    tt = sum( DoubleDouble.Make( s1, s2 ), 2 );
                    [ t1( jdxs ), t2( jdxs ) ] = DoubleDouble.DDPlusDD( v.v1( jdxs, j ), v.v2( jdxs, j ), -tt.v1, -tt.v2 );
                    [ L.v1( jdxs, j ), L.v2( jdxs, j ) ] = DoubleDouble.DDDividedByDD( t1( jdxs ), t2( jdxs ), x1( j ), x2( j ) );
                end
            end
            if nargin < 2 || ~strcmp( type, 'vector_d' )
                d = diag( d );
            else
                d = d.';
            end
        end
        
        function [ v, d ] = eig( x )
            [ v, d ] = eig( x.v1 );
            v = DoubleDouble( v );
            d = DoubleDouble( diag( d ) );
            C = length( d );
            I = eye( C, 'DoubleDouble' );
            for c = 1 : C
                vi = DoubleDouble.Make( v.v1( :, c ), v.v2( :, c ) );
                dii = DoubleDouble.Make( d.v1( c, 1 ), d.v2( c, 1 ) );
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
                d.v1( c, 1 ) = dii.v1;
                d.v2( c, 1 ) = dii.v2;
                v.v1( :, c ) = vi.v1;
                v.v2( :, c ) = vi.v2;
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

            w = DoubleDouble.zeros( K, 1 );

            for k = 1 : K

                j = max( 1, k + 1 - N ) : min( k, M );
                i = k - j + 1;

                wk = DoubleDouble.Dot( DoubleDouble.Make( u.v1( j ), u.v2( j ) ), DoubleDouble.Make( v.v1( i ), v.v2( i ) ) );
                w.v1( k ) = wk.v1;
                w.v2( k ) = wk.v2;

            end

            if RowVector
                w = w.';
            end

        end

    end

    methods ( Static )
        function v = IsEqualWithExpansion( a, b, varargin )
            v = a == b;
            v = all( v(:) );
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
            v = DoubleDouble.Make( ones( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end
        
        function v = zeros( varargin )
            v = DoubleDouble.Make( zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end
        
        function v = eye( varargin )
            v = DoubleDouble.Make( eye( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end
        
        function v = nan( varargin )
            v = DoubleDouble.Make( nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ) );
        end
        
        function v = inf( varargin )
            v = DoubleDouble.Make( inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ) );
        end
        
        function v = rand( varargin )
            t = rand( varargin{:}, 'double' );
            v = DoubleDouble.Make( t, eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ) );
        end
        
        function v = randn( varargin )
            t = randn( varargin{:}, 'double' );
            v = DoubleDouble.Make( t, eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ) );
        end
        
        function v = randi( imax, varargin )
            v = DoubleDouble.Make( randi( imax, varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end
        
        function v = Plus( a, b )
            if isa( a, 'DoubleDouble' )
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDPlusDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DDPlusDouble( a.v1, a.v2, double( b ) );
                end
            else
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDPlusDouble( b.v1, b.v2, double( a ) );
                else
                    [ x1, x2 ] = DoubleDouble.DoublePlusDouble( double( a ), double( b ) );
                end
            end
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = Minus( a, b )
            if isa( a, 'DoubleDouble' )
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDPlusDD( a.v1, a.v2, -b.v1, -b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DDPlusDouble( a.v1, a.v2, -double( b ) );
                end
            else
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDPlusDouble( -b.v1, -b.v2, double( a ) );
                else
                    [ x1, x2 ] = DoubleDouble.DoublePlusDouble( double( a ), -double( b ) );
                end
            end
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = Times( a, b )
            if isa( a, 'DoubleDouble' )
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDTimesDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DDTimesDouble( a.v1, a.v2, double( b ) );
                end
            else
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDTimesDouble( b.v1, b.v2, double( a ) );
                else
                    [ x1, x2 ] = DoubleDouble.DoubleTimesDouble( double( a ), double( b ) );
                end
            end
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = MTimes( a, b )
            [ R, c ] = size( a );
            [ r, C ] = size( b );
            if ( ( R == 1 ) && ( c == 1 ) ) || ( ( r == 1 ) && ( C == 1 ) )
                v = DoubleDouble.Times( a, b );
                return
            end
            v = DoubleDouble.Make( zeros( R, C ), zeros( R, C ) );
            if isa( b, 'DoubleDouble' )
                for c = 1 : C
                    t = DoubleDouble.Sum( a .* DoubleDouble.Make( b.v1( :, c ).', b.v2( :, c ).' ), 2 );
                    v.v1( :, c ) = t.v1;
                    v.v2( :, c ) = t.v2;
                end
            else
                for c = 1 : C
                    t = DoubleDouble.Sum( a .* b( :, c ).', 2 );
                    v.v1( :, c ) = t.v1;
                    v.v2( :, c ) = t.v2;
                end
            end
        end
        
        function v = RDivide( a, b )
            if isa( a, 'DoubleDouble' )
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDDividedByDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DDDividedByDouble( a.v1, a.v2, double( b ) );
                end
            else
                if isa( b, 'DoubleDouble' )
                    da = double( a );
                    [ x1, x2 ] = DoubleDouble.DDDividedByDD( da, zeros( size( da ) ), b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DoubleDividedByDouble( double( a ), double( b ) );
                end
            end
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = LDivide( b, a )
            if isa( a, 'DoubleDouble' )
                if isa( b, 'DoubleDouble' )
                    [ x1, x2 ] = DoubleDouble.DDDividedByDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DDDividedByDouble( a.v1, a.v2, double( b ) );
                end
            else
                if isa( b, 'DoubleDouble' )
                    da = double( a );
                    [ x1, x2 ] = DoubleDouble.DDDividedByDD( da, zeros( size( da ) ), b.v1, b.v2 );
                else
                    [ x1, x2 ] = DoubleDouble.DoubleDividedByDouble( double( a ), double( b ) );
                end
            end
            v = DoubleDouble.Make( x1, x2 );
        end
        
        function v = MLDivide( a, v )
            [ Ra, Ca ] = size( a );
            [ Rv, Cv ] = size( v );
            if ( ( Ra == 1 ) && ( Ca == 1 ) ) || ( ( Rv == 1 ) && ( Cv == 1 ) )
                v = DoubleDouble.LDivide( a, v );
                return
            end
            if ~isa( v, 'DoubleDouble' )
                v = DoubleDouble( v );
            end
            assert( Ra == Rv );
            if Ra ~= Ca
                [ q, a ] = qr( a );
                v = q' * v;
                v = BackSubstitution( v, a );
                return
            end
            if DoubleDouble.IsEqualWithExpansion( triu( a, 1 ), 0 )
                % Lower triangular
                v = ForwardElimination( v, a );
                return
            elseif DoubleDouble.IsEqualWithExpansion( tril( a, -1 ), 0 )
                % Upper triangular
                v = BackSubstitution( v, a );
                return
            elseif DoubleDouble.IsEqualWithExpansion( a, a' )
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
            v.v1 = v.v1( p, : );
            v.v2 = v.v2( p, : );
            v = ForwardElimination( v, L );

            % Back substitution
            v = BackSubstitution( v, U );
        end
        
        function v = MRDivide( v, a )
            v = DoubleDouble.MLDivide( a', v' )';
        end
        
        function [ v, Indices ] = Sort( v, dim, cm, varargin )
            if nargin < 2 || isempty( dim )
                dim = find( size( v.v1 ) > 1, 1 );
                if isempty( dim )
                    dim = 1;
                end
            end
            if nargin < 3 || isempty( cm )
                cm = 'auto';
            end
            if isa( v, 'DoubleDouble' )
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
                Blocks{ dim } = Size( dim );
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
                v       = DoubleDouble.Make( cell2mat( xv1 ), cell2mat( xv2 ) );
            else
                if nargout > 1
                    [ v, Indices ] = sort( v, dim, 'ComparisonMethod', cm, varargin{:} );
                else
                    v = sort( v, dim, 'ComparisonMethod', cm, varargin{:} );
                end
                v = DoubleDouble( v );
            end
        end
        
        function s = Sum( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                end
            else
                if nargin < 2 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, x{ i } );
                end
            end
        end
        
        function c = CumSum( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 2 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end
        
        function c = Diff( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = DoubleDouble.Make( x1{ i }, x2{ i } );
                    d = DoubleDouble.Minus( t, s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            else
                if nargin < 2 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = [];
                c2{1} = [];
                for i = 2 : Length
                    t = x{ i };
                    d = DoubleDouble.Minus( t, s );
                    c1{i} = d.v1;
                    c2{i} = d.v2;
                    s = t;
                end
            end
            c = DoubleDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end
        
        function s = Prod( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = DoubleDouble.Make( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = DoubleDouble.Times( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                end
            else
                if nargin < 2 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    Size = max( 1, Size );
                    s = DoubleDouble.Make( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = DoubleDouble.Times( s, x{ i } );
                end
            end
        end
        
        function c = CumProd( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Times( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 2 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = DoubleDouble.Times( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
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
                    v = DoubleDouble;
                end
            end
        end
        
        function [ s, i ] = Max( a, b, dim )
            if isempty( b )
                if isempty( a )
                    s = DoubleDouble;
                    i = [];
                    return
                end
                if isa( a, 'DoubleDouble' )
                    if nargin < 3 || isempty( dim )
                        dim = find( size( a.v1 ) > 1, 1 );
                        if isempty( dim )
                            dim = 1;
                        end
                    end
                    Size = size( a.v1 );
                    Length = Size( dim );
                    Blocks = num2cell( Size );
                    Blocks{ dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1, Blocks{:} );
                    x2 = mat2cell( a.v2, Blocks{:} );
                    s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                    Size( dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = DoubleDouble.Max( DoubleDouble.Make( x1{ j }, x2{ j } ), s );
                        i( ii ) = j;
                    end
                else
                    if nargin < 3 || isempty( dim )
                        dim = find( size( a ) > 1, 1 );
                        if isempty( dim )
                            dim = 1;
                        end
                    end
                    Size = size( a );
                    Length = Size( dim );
                    Blocks = num2cell( Size );
                    Blocks{ dim } = ones( Length, 1 );
                    x = mat2cell( a, Blocks{:} );
                    s = x{ 1 };
                    for j = 2 : Length
                        s = DoubleDouble.Max( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'DoubleDouble' )
                    a = DoubleDouble( a );
                end
                if ~isa( b, 'DoubleDouble' )
                    b = DoubleDouble( b );
                end
                [ a, b ] = DoubleDouble.ExpandSingleton( a, b );
                i = ( a.v1 > b.v1 ) | ( ( a.v1 == b.v1 ) & ( a.v2 > b.v2 ) );
                s = b;
                s.v1( i ) = a.v1( i );
                s.v2( i ) = a.v2( i );
            end
        end
        
        function c = CumMax( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 3 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Max( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 3 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = DoubleDouble.Max( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end
        
        function [ s, i ] = Min( a, b, dim )
            if isempty( b )
                if isa( a, 'DoubleDouble' )
                    if nargin < 3 || isempty( dim )
                        dim = find( size( a.v1 ) > 1, 1 );
                        if isempty( dim )
                            dim = 1;
                        end
                    end
                    Size = size( a.v1 );
                    Length = Size( dim );
                    Blocks = num2cell( Size );
                    Blocks{ dim } = ones( Length, 1 );
                    x1 = mat2cell( a.v1, Blocks{:} );
                    x2 = mat2cell( a.v2, Blocks{:} );
                    s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                    Size( dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = DoubleDouble.Min( DoubleDouble.Make( x1{ j }, x2{ j } ), s );
                        i( ii ) = j;
                    end
                else
                    if nargin < 3 || isempty( dim )
                        dim = find( size( a ) > 1, 1 );
                        if isempty( dim )
                            dim = 1;
                        end
                    end
                    Size = size( a );
                    Length = Size( dim );
                    Blocks = num2cell( Size );
                    Blocks{ dim } = ones( Length, 1 );
                    x = mat2cell( a, Blocks{:} );
                    s = x{ 1 };
                    for j = 2 : Length
                        s = DoubleDouble.Min( s, x{ j } );
                    end
                end
            else
                if ~isa( a, 'DoubleDouble' )
                    a = DoubleDouble( a );
                end
                if ~isa( b, 'DoubleDouble' )
                    b = DoubleDouble( b );
                end
                [ a, b ] = DoubleDouble.ExpandSingleton( a, b );
                i = ( a.v1 < b.v1 ) | ( ( a.v1 == b.v1 ) & ( a.v2 < b.v2 ) );
                s = b;
                s.v1( i ) = a.v1( i );
                s.v2( i ) = a.v2( i );
            end
        end
        
        function c = CumMin( v, dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 3 || isempty( dim )
                    dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.Make( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Min( s, DoubleDouble.Make( x1{ i }, x2{ i } ) );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            else
                if nargin < 3 || isempty( dim )
                    dim = find( size( v ) > 1, 1 );
                    if isempty( dim )
                        dim = 1;
                    end
                end
                Size = size( v );
                Length = Size( dim );
                if Length == 0
                    c = DoubleDouble.Make( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                c1 = cell( size( x ) );
                c2 = cell( size( x ) );
                c1{1} = s;
                c2{1} = zeros( size( s ) );
                for i = 2 : Length
                    s = DoubleDouble.Min( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.Make( cell2mat( c1 ), cell2mat( c2 ) );
        end
                
        function v = Dot( a, b, dim )
            if nargin < 3
                dim = [];
            end
            if ( length( a ) == numel( a ) ) && ( length( b ) == numel( b ) )
                a = Vec( a );
                b = Vec( b );
            end
            v = DoubleDouble.Sum( DoubleDouble.Times( a, b ), dim );
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
            v.v1 = v.v1(:);
            v.v2 = v.v2(:);
        end

        function v = TimesPowerOf2( v, b )
            assert( isa( b, 'double' ) );
            v.v1 = v.v1 .* b;
            v.v2 = v.v2 .* b;
        end

        function v = ForwardElimination( v, L )
            % For lower triangular L, x = ForwardElimination( b, L ) solves L*x = b.
            [ m, n ] = size( L );
            mn = min( m, n );
            [ vm, vn ] = size( v );
            if vm < n
                v = [ v; zeros( n - vm, vn ) ];
            elseif vm > n
                v.v1( (n+1):vm, : ) = [];
                v.v2( (n+1):vm, : ) = [];
            end
            if isa( L, 'DoubleDouble' )
                [ v.v1( 1, : ), v.v2( 1, : ) ] = DoubleDouble.DDDividedByDD( v.v1( 1, : ), v.v2( 1, : ), L.v1( 1, 1 ), L.v2( 1, 1 ), true );
                for k = 2 : mn
                   j = 1 : k - 1;
                   [ t1, t2 ] = DoubleDouble.DDTimesDD( v.v1( j, : ), v.v2( j, : ), L.v1( k, j ).', L.v2( k, j ).' );
                   t = DoubleDouble.Sum( DoubleDouble.Make( t1, t2 ), 1 );
                   [ t1, t2 ] = DoubleDouble.DDPlusDD( v.v1( k, : ), v.v2( k, : ), -t.v1, -t.v2 );
                   [ v.v1( k, : ), v.v2( k, : ) ] = DoubleDouble.DDDividedByDD( t1, t2, L.v1( k, k ), L.v2( k, k ), true );
                end
            else
                [ v.v1( 1, : ), v.v2( 1, : ) ] = DoubleDouble.DDDividedByDouble( v.v1( 1, : ), v.v2( 1, : ), L( 1, 1 ), true );
                for k = 2 : mn
                   j = 1 : k - 1;
                   [ t1, t2 ] = DoubleDouble.DDTimesDouble( v.v1( j, : ), v.v2( j, : ), L( k, j ).' );
                   t = DoubleDouble.Sum( DoubleDouble.Make( t1, t2 ), 1 );
                   [ t1, t2 ] = DoubleDouble.DDPlusDD( v.v1( k, : ), v.v2( k, : ), -t.v1, -t.v2 );
                   [ v.v1( k, : ), v.v2( k, : ) ] = DoubleDouble.DDDividedByDouble( t1, t2, L( k, k ), true );
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
                v.v1( (n+1):vm, : ) = [];
                v.v2( (n+1):vm, : ) = [];
            end
            if isa( U, 'DoubleDouble' )
                [ v.v1( mn, : ), v.v2( mn, : ) ] = DoubleDouble.DDDividedByDD( v.v1( mn, : ), v.v2( mn, : ), U.v1( mn, mn ), U.v2( mn, mn ), true );
                for k = mn - 1 : -1 : 1
                   j = k + 1 : n;
                   [ t1, t2 ] = DoubleDouble.DDTimesDD( v.v1( j, : ), v.v2( j, : ), U.v1( k, j ).', U.v2( k, j ).' );
                   t = DoubleDouble.Sum( DoubleDouble.Make( t1, t2 ), 1 );
                   [ t1, t2 ] = DoubleDouble.DDPlusDD( v.v1( k, : ), v.v2( k, : ), -t.v1, -t.v2 );
                   [ v.v1( k, : ), v.v2( k, : ) ] = DoubleDouble.DDDividedByDD( t1, t2, U.v1( k, k ), U.v2( k, k ), true );
                end
            else
                [ v.v1( mn, : ), v.v2( mn, : ) ] = DoubleDouble.DDDividedByDouble( v.v1( mn, : ), v.v2( mn, : ), U( mn, mn ), true );
                for k = mn - 1 : -1 : 1
                   j = k + 1 : n;
                   [ t1, t2 ] = DoubleDouble.DDTimesDouble( v.v1( j, : ), v.v2( j, : ), U( k, j ).' );
                   t = DoubleDouble.Sum( DoubleDouble.Make( t1, t2 ), 1 );
                   [ t1, t2 ] = DoubleDouble.DDPlusDD( v.v1( k, : ), v.v2( k, : ), -t.v1, -t.v2 );
                   [ v.v1( k, : ), v.v2( k, : ) ] = DoubleDouble.DDDividedByDouble( t1, t2, U( k, k ), true );
                end
            end
        end
        
        function v = SinTaylor( v )
            Threshhold = 0.5 .* abs( v.v1(:) ) .* DoubleDouble.eps.v1;
            x = - v .* v;
            r = v;
            for i = 1 : 2 : DoubleDouble.NInverseFactorial
                r = r .* x;
                t = r .* DoubleDouble.Make( DoubleDouble.InverseFactorial( i, 1 ), DoubleDouble.InverseFactorial( i, 2 ) );
                v = v + t;
                if all( abs( t.v1(:) ) <= Threshhold )
                    break
                end
            end
        end
        
        function [ sin_v, cos_v ] = SinCosTaylor( v )
            sin_v = SinTaylor( v );
            cos_v = sqrt( 1 - sin_v .* sin_v );
        end
    end
    
    methods ( Static, Access = private )
        function v = Make( a1, a2 )
            v = DoubleDouble;
            v.v1 = a1;
            v.v2 = a2;
        end
        
        function [ s1, s2 ] = Normalize( a1, a2 )
            s1 = a1 + a2;
            t = s1 - a1;
            s2 = a2 - t;
        end
        
        function [ s1, s2 ] = DDPlusDD( a1, a2, b1, b2 )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = DoubleDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            [ s1, s2 ] = DoubleDouble.DoublePlusDouble( a1, b1 );
            [ t1, t2 ] = DoubleDouble.DoublePlusDouble( a2, b2 );
            s2 = s2 + t1;
            [ s1, s2 ] = DoubleDouble.Normalize( s1, s2 );
            s2 = s2 + t2;
            [ s1, s2 ] = DoubleDouble.Normalize( s1, s2 );
        end

        function [ s1, s2 ] = DDPlusDouble( a1, a2, b )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = DoubleDouble.ExpandSingleton( a1, a2, b );
            end
            [ s1, s2 ] = DoubleDouble.DoublePlusDouble( a1, b );
            s2 = s2 + a2;
            [ s1, s2 ] = DoubleDouble.Normalize( s1, s2 );
        end

        function [ s1, s2 ] = DoublePlusDouble( a, b )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a, b ] = DoubleDouble.ExpandSingleton( a, b );
            end
            s1 = a + b;
            bb = s1 - a;
            t11 = s1 - bb;
            t2 = b - bb;
            t1 = a - t11;
            s2 = t1 + t2;
        end

        function [ p1, p2 ] = DDTimesDD( a1, a2, b1, b2 )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = DoubleDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            [ p1, p2 ] = DoubleDouble.DoubleTimesDouble( a1, b1 );
            t = a1 .* b2 + a2 .* b1;
            p2 = p2 + t;
            [ p1, p2 ] = DoubleDouble.Normalize( p1, p2 );
        end

        function [ p1, p2 ] = DDTimesDouble( a1, a2, b )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = DoubleDouble.ExpandSingleton( a1, a2, b );
            end
            [ p1, p2 ] = DoubleDouble.DoubleTimesDouble( a1, b );
            p2 = p2 + a2 .* b;
            [ p1, p2 ] = DoubleDouble.Normalize( p1, p2 );
        end

        function [ p1, p2 ] = DoubleTimesDouble( a, b )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a, b ] = DoubleDouble.ExpandSingleton( a, b );
            end
            p1 = a .* b;
            [ a1, a2 ] = DoubleDouble.Split( a );
            [ b1, b2 ] = DoubleDouble.Split( b );
            t1 = a1 .* b1 - p1;
            t2 = t1 + a1 .* b2 + a2 .* b1;
            p2 = t2 + a2 .* b2;
        end
        
        function [ r1, r2 ] = DDDividedByDD( a1, a2, b1, b2, AnySolutionWillDo )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b1, b2 ] = DoubleDouble.ExpandSingleton( a1, a2, b1, b2 );
            end
            if nargin < 5
                AnySolutionWillDo = false;
            end
            q1 = a1 ./ b1;
            [ p1, p2 ] = DoubleDouble.DDTimesDouble( b1, b2, q1 );
            [ r1, r2 ] = DoubleDouble.DDPlusDD( a1, a2, -p1, -p2 );
            q2 = r1 ./ b1;
            [ p1, p2 ] = DoubleDouble.DDTimesDouble( b1, b2, q2 );
            [ r1, ~  ] = DoubleDouble.DDPlusDD( r1, r2, -p1, -p2 );
            q3 = r1 ./ b1;
            [ q1, q2 ] = DoubleDouble.Normalize( q1, q2 );
            [ r1, r2 ] = DoubleDouble.DDPlusDD( q1, q2, q3, zeros( size( q3 ) ) );
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
                    a1Select( a1SelectSelect ) = 0;
                end
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
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
                a1Select( a1SelectSelect ) = NaN;
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
            end
        end
        
        function [ r1, r2 ] = DDDividedByDouble( a1, a2, b, AnySolutionWillDo )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a1, a2, b ] = DoubleDouble.ExpandSingleton( a1, a2, b );
            end
            if nargin < 4
                AnySolutionWillDo = false;
            end
            r1 = a1 ./ b;
            [ p1, p2 ] = DoubleDouble.DoubleTimesDouble( r1, b );
            [ s, e ] = DoubleDouble.DoublePlusDouble( a1, -p1 );
            e = e + a2;
            e = e - p2;
            t = s + e;
            r2 = t ./ b;
            [ r1, r2 ] = DoubleDouble.Normalize( r1, r2 );
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
                    a1Select( a1SelectSelect ) = 0;
                end
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
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
                a1Select( a1SelectSelect ) = NaN;
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
            end
        end
        
        function [ r1, r2 ] = DoubleDividedByDouble( a, b, AnySolutionWillDo )
            if DoubleDouble.SingletonExpansionNotSupported
                [ a, b ] = DoubleDouble.ExpandSingleton( a, b );
            end
            if nargin < 3
                AnySolutionWillDo = false;
            end
            r1 = a ./ b;
            [ p1, p2 ] = DoubleDouble.DoubleTimesDouble( r1, b );
            [ s, e ] = DoubleDouble.DoublePlusDouble( a, -p1 );
            e = e - p2;
            t = s + e;
            r2 = t ./ b;
            [ r1, r2 ] = DoubleDouble.Normalize( r1, r2 );
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
                    a1Select( a1SelectSelect ) = 0;
                end
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
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
                a1Select( a1SelectSelect ) = NaN;
                r1( Select ) = a1Select;
                r2( Select ) = a1Select;
            end
        end
        
        function [ a1, a2 ] = Split( a )
            if isreal( a )
                Select = ( a > 6.69692879491417e+299 ) | ( a < -6.69692879491417e+299 ); % 2^996
                a( Select ) = a( Select ) * 3.7252902984619140625e-09; % 2^(-28)
                t1 = 134217729.0 * a; % 2^27 + 1
                t2 = t1 - a;
                a1 = t1 - t2;
                a2 = a - a1;
                a1( Select ) = a1( Select ) * 268435456.0; % 2^28
                a2( Select ) = a2( Select ) * 268435456.0; % 2^28
            else
                [ r1, r2 ] = DoubleDouble.Split( real( a ) );
                [ i1, i2 ] = DoubleDouble.Split( imag( a ) );
                a1 = complex( r1, i1 );
                a2 = complex( r2, i2 );
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
