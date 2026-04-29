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

        expRescale = 9;
        logSteps = 2;

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
            v = ToRand( DoubleDouble.zeros( varargin{:} ) );
        end

        function v = randn( varargin )
            v = ToRandn( DoubleDouble.zeros( varargin{:} ) );
        end

        function v = randi( imax, varargin )
            v = DoubleDouble.MakeStatic( randi( imax, varargin{:}, 'double' ), zeros( varargin{:}, 'double' ) );
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
