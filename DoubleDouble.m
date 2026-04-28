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

classdef DoubleDouble < BaseExtDouble

    properties ( Constant, GetAccess = public )
        zero = DoubleDouble.MakeStatic( 0, 0 );
        one = DoubleDouble.MakeStatic( 1, 0 );
        eps = DoubleDouble.MakeStatic( 4.93038065763132e-32, 0 );
        pi = DoubleDouble.MakeStatic( 3.141592653589793116e+00, 1.224646799147353207e-16 );
    end

    properties ( Constant, GetAccess = public )

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

        piT2   = DoubleDouble.MakeStatic( 6.283185307179586232e+00,  2.449293598294706414e-16 );
        piD2   = DoubleDouble.MakeStatic( 1.570796326794896558e+00,  6.123233995736766036e-17 );
        piD16  = DoubleDouble.MakeStatic( 1.963495408493620697e-01,  7.654042494670957545e-18 );

        log_2  = DoubleDouble.MakeStatic( 6.931471805599452862e-01,  2.319046813846299558e-17 );
        log_10 = DoubleDouble.MakeStatic( 2.302585092994045901e+00, -2.170756223382249351e-16 );

        SinTable = [
            0, 0;
            1.950903220161282758e-01, -7.991079068461731263e-18;
            3.826834323650897818e-01, -1.005077269646158761e-17;
            5.555702330196021776e-01,  4.709410940561676821e-17;
            7.071067811865475727e-01, -4.833646656726456726e-17;
            ];

        CosTable = [
            1, 0;
            9.807852804032304306e-01,  1.854693999782500573e-17;
            9.238795325112867385e-01,  1.764504708433667706e-17;
            8.314696123025452357e-01,  1.407385698472802389e-18;
            7.071067811865475727e-01, -4.833646656726456726e-17;
            ];
    end

    methods

        function v = Promote( ~, a )
            v = DoubleDouble( a );
        end

        function v = Make( ~, a1, a2 )
            v = DoubleDouble.MakeStatic( a1, a2 );
        end

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
            elseif isa( in, 'BaseExtDouble' )
                [ v.v1, v.v2 ] = ToSumOfDoubles( in );
            else
                v.v1 = double( in );
                v.v2 = zeros( size( in ) );
            end
        end

    end

    methods ( Static )

        %% I am happy that these should not be in the base class.

        function v = ones( varargin )
            v = DoubleDouble.MakeStatic( ones( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = zeros( varargin )
            v = DoubleDouble.MakeStatic( zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = eye( varargin )
            v = DoubleDouble.MakeStatic( eye( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = nan( varargin )
            v = DoubleDouble.MakeStatic( nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ) );
        end

        function v = inf( varargin )
            v = DoubleDouble.MakeStatic( inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ) );
        end

        function v = rand( varargin )
            t = rand( varargin{:}, 'double' );
            v = DoubleDouble.MakeStatic( t, eps( t ) .* ( rand( varargin{:}, 'double' ) - 0.5 ) );
        end

        function v = randi( imax, varargin )
            v = DoubleDouble.MakeStatic( randi( imax, varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end



        %% The following should have the implementations in the base class (non-static, sealed), and the static functions should be trivial wrappers.

        function v = randn( varargin )
            Size = [ varargin{ : } ];
            if isempty( Size )
                Size = [ 1, 1 ];
            elseif isscalar( Size )
                Size = [ Size, Size ];
            end
            N = prod( Size );
            M = 2 * ceil( 0.5 * N );
            U = DoubleDouble.rand( M / 2, 1 );
            V = DoubleDouble.rand( M / 2, 1 );
            R = sqrt( -2 * log( U ) );
            Theta = 2 * DoubleDouble.pi * V;
            [ S, C ] = sincos( Theta );
            Z1 = R .* C;
            Z2 = R .* S;
            v = reshape( [ Z1, Z2 ].', [ M, 1 ] );
            v.v1 = v.v1( 1 : N );
            v.v2 = v.v2( 1 : N );
            v = reshape( v, Size );
        end









        %% The following should be non-static sealed functions in the base class.

        function s = Sum( v, Dim )
            if isa( v, 'DoubleDouble' )
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
                    s = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    s = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, x{ i } );
                end
            end
        end

        function c = CumSum( v, Dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Plus( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
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
                    s = DoubleDouble.Plus( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.MakeStatic( cell2mat( c1 ), cell2mat( c2 ) );
        end


        function s = Prod( v, Dim )
            if isa( v, 'DoubleDouble' )
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
                    s = DoubleDouble.MakeStatic( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                for i = 2 : Length
                    s = DoubleDouble.Times( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    s = DoubleDouble.MakeStatic( ones( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x = mat2cell( v, Blocks{:} );
                s = x{ 1 };
                for i = 2 : Length
                    s = DoubleDouble.Times( s, x{ i } );
                end
            end
        end

        function c = CumProd( v, Dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 2 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Times( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
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
                    s = DoubleDouble.Times( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.MakeStatic( cell2mat( c1 ), cell2mat( c2 ) );
        end


        function [ s, i ] = Max( a, b, Dim )
            if isempty( b )
                if isempty( a )
                    s = DoubleDouble;
                    i = [];
                    return
                end
                if isa( a, 'DoubleDouble' )
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
                    s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = DoubleDouble.Max( DoubleDouble.MakeStatic( x1{ j }, x2{ j } ), s );
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

        function c = CumMax( v, Dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Max( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
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
                    s = DoubleDouble.Max( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.MakeStatic( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function [ s, i ] = Min( a, b, Dim )
            if isempty( b )
                if isa( a, 'DoubleDouble' )
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
                    s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                    Size( Dim ) = 1;
                    i = ones( Size );
                    for j = 2 : Length
                        [ s, ii ] = DoubleDouble.Min( DoubleDouble.MakeStatic( x1{ j }, x2{ j } ), s );
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

        function c = CumMin( v, Dim )
            if isa( v, 'DoubleDouble' )
                if nargin < 3 || isempty( Dim )
                    Dim = find( size( v.v1 ) > 1, 1 );
                    if isempty( Dim )
                        Dim = 1;
                    end
                end
                Size = size( v.v1 );
                Length = Size( Dim );
                if Length == 0
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
                    return
                end
                Blocks = num2cell( Size );
                Blocks{ Dim } = ones( Length, 1 );
                x1 = mat2cell( v.v1, Blocks{:} );
                x2 = mat2cell( v.v2, Blocks{:} );
                s = DoubleDouble.MakeStatic( x1{ 1 }, x2{ 1 } );
                c1 = cell( size( x1 ) );
                c2 = cell( size( x2 ) );
                c1{1} = s.v1;
                c2{1} = s.v2;
                for i = 2 : Length
                    s = DoubleDouble.Min( s, DoubleDouble.MakeStatic( x1{ i }, x2{ i } ) );
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
                    c = DoubleDouble.MakeStatic( zeros( Size ), zeros( Size ) );
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
                    s = DoubleDouble.Min( s, x{ i } );
                    c1{i} = s.v1;
                    c2{i} = s.v2;
                end
            end
            c = DoubleDouble.MakeStatic( cell2mat( c1 ), cell2mat( c2 ) );
        end

        function v = Dot( a, b, Dim )
            if nargin < 3
                Dim = [];
            end
            if ( length( a ) == numel( a ) ) && ( length( b ) == numel( b ) )
                a = Vec( a );
                b = Vec( b );
            end
            v = DoubleDouble.Sum( DoubleDouble.Times( a, b ), Dim );
        end

    end


    methods ( Static, Access = public )

        function v = MakeStatic( a1, a2 )
            v = DoubleDouble;
            v.v1 = a1;
            v.v2 = a2;
        end








    end
end
