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
                C = cell( 1, 4 );
                [ C{ : } ] = ToSumOfDoubles( in );
                if length( C ) < 2 || isempty( C{ 2 } ); C{ 2 } = zeros( size( C{ 1 } ) ); end
                if length( C ) < 3 || isempty( C{ 3 } ); C{ 3 } = zeros( size( C{ 1 } ) ); end
                if length( C ) < 4 || isempty( C{ 4 } ); C{ 4 } = zeros( size( C{ 1 } ) ); end
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );

            else
                v.v1 = DoubleDouble.MakeStatic( double( in ), zeros( size( in ) ) );
                v.v2 = DoubleDouble.MakeStatic( zeros( size( in ) ), zeros( size( in ) ) );
            end
        end

        function v = Promote( ~, v )
            v = QuadDouble( v );
        end

        function v = Make( ~, a1, a2 )
            v = QuadDouble;
            if ~isa( a1, 'DoubleDouble' ) && ~isempty( a1 )
                a1 = DoubleDouble( a1 );
            end
            if ~isa( a2, 'DoubleDouble' ) && ~isempty( a2 )
                a2 = DoubleDouble( a2 );
            end
            v.v1 = a1;
            v.v2 = a2;
        end

    end

    methods ( Static )

        function v = ones( varargin )
            v = QuadDouble.MakeStatic( ones( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = zeros( varargin )
            v = QuadDouble.MakeStatic( zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = eye( varargin )
            v = QuadDouble.MakeStatic( eye( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

        function v = nan( varargin )
            v = QuadDouble.MakeStatic( nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ), nan( varargin{:}, 'double' ) );
        end

        function v = inf( varargin )
            v = QuadDouble.MakeStatic( inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ), inf( varargin{:}, 'double' ) );
        end

        function v = randi( imax, varargin )
            v = QuadDouble.MakeStatic( randi( imax, varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
        end

    end

    methods ( Static, Access = public )

        function v = MakeStatic( a1, a2, a3, a4 )
            v = QuadDouble;
            if nargin < 2
                a2 = zeros( size( a1 ) );
            end
            v.v1 = DoubleDouble.MakeStatic( a1, a2 );
            if nargin < 3
                a3 = zeros( size( a1 ) );
            end
            if nargin < 4
                a4 = zeros( size( a1 ) );
            end
            v.v2 = DoubleDouble.MakeStatic( a3, a4 );
        end

    end

    methods ( Access = protected )

        function v = Plus( a, b )
            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QuadDouble.QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( s0, s1, s2, s3 );
        end

        function v = Times( a, b )
            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QuadDouble.QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( s0, s1, s2, s3 );
        end

        function v = RDivide( a, b )
            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QuadDouble.QDDivQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( s0, s1, s2, s3 );
        end

    end

    methods ( Static, Access = private )

        function [ s1, s2 ] = Normalize( a1, a2 )
            s1 = a1 + a2;
            t = s1 - a1;
            s2 = a2 - t;
        end

        function v = PromoteStatic( a )
            v = QuadDouble( a );
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


        function [ p1, p2 ] = TwoProd( a, b )
            if QuadDouble.SingletonExpansionNotSupported
                [ a, b ] = QuadDouble.ExpandSingleton( a, b );
            end
            p1 = a .* b;
            [ a1, a2 ] = Split( a );
            [ b1, b2 ] = Split( b );
            p2 = ( ( a1 .* b1 - p1 ) + a1 .* b2 + a2 .* b1 ) + a2 .* b2;
        end

        function [ p1, p2 ] = TwoSqr( a )
            p1 = a .* a;
            [ a1, a2 ] = Split( a );
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
    end

    methods ( Static, Access = public )
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

    end

end
