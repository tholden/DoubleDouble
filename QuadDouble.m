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
                v.v1 = DoubleDouble.empty;
                v.v2 = DoubleDouble.empty;
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = QuadDouble.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDouble' ) || isa( in, 'QuadDoubleSlow' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            elseif isa( in, 'DoubleDouble' )
                v.v1 = in;
                v.v2 = DoubleDouble.zeros( size( in ) );
            elseif isa( in, 'BaseExtDouble' )
                C = cell( 1, 4 );
                [ C{ : } ] = ToSumOfDoubles( in );
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );
            else
                v.v1 = DoubleDouble( in );
                v.v2 = DoubleDouble.zeros( size( in ) );
            end
        end

        function v = Promote( ~, v )
            v = QuadDouble( v );
        end

        function n = PromotionOrder( ~ )
            n = 2.5;
        end

    end

    methods ( Static )

        function v = ones( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.ones( varargin{:} ), DoubleDouble.zeros( varargin{:} ) );
        end

        function v = zeros( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.zeros( varargin{:} ), DoubleDouble.zeros( varargin{:} ) );
        end

        function v = eye( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.eye( varargin{:} ), DoubleDouble.zeros( varargin{:} ) );
        end

        function v = NaN( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.NaN( varargin{:} ), DoubleDouble.NaN( varargin{:} ) );
        end

        function v = Inf( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.Inf( varargin{:} ), DoubleDouble.Inf( varargin{:} ) );
        end

        function v = rand( varargin )
            v = ToRand( QuadDouble.zeros( varargin{:} ) );
        end

        function v = randn( varargin )
            v = ToRandn( QuadDouble.zeros( varargin{:} ) );
        end

        function v = randi( imax, varargin )
            v = QuadDouble.MakeStatic( DoubleDouble( randi( imax, varargin{:}, 'double' ) ), DoubleDouble.zeros( varargin{:} ) );
        end

    end

    methods ( Static, Access = ?BaseExtDouble )

        function v = MakeStatic( a1, a2 )
            v = QuadDouble;
            if ~isa( a1, 'DoubleDouble' )
                a1 = DoubleDouble( a1 );
            end
            if ~isa( a2, 'DoubleDouble' )
                a2 = DoubleDouble( a2 );
            end
            v.v1 = a1;
            v.v2 = a2;
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDouble.MakeStatic( a1, a2 );
        end

        function v = Plus( a, b )
            % [ a, b ] = BaseExtDouble.JointPromotion( a, b );
            % if PromotionOrder( a ) == PromotionOrder( b )
            %     [ x1, x2 ] = QDPlusQD( a.v1, a.v2, b.v1, b.v2 );
            %     v = a.Make( x1, x2 );
            % elseif PromotionOrder( a ) > PromotionOrder( b )
            %     [ x1, x2 ] = QuadDouble.QDPlusDD( a.v1, a.v2, Promote( a.v1, b ) );
            %     v = a.Make( x1, x2 );
            % else % PromotionOrder( b ) > PromotionOrder( a )
            %     [ x1, x2 ] = QuadDouble.QDPlusDD( b.v1, b.v2, Promote( b.v1, a ) );
            %     v = b.Make( x1, x2 );
            % end

            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

        function v = Times( a, b )
            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

        function v = RDivide( a, b )
            a = QuadDouble.PromoteStatic( a );
            b = QuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDDividedByQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

    end

    methods ( Static, Access = private )

        function v = PromoteStatic( a ) % TODO Remove
            v = QuadDouble( a );
        end

    end

end
