classdef DoubleDouble < BaseDoubleDouble & ExtDouble

    properties ( Constant, GetAccess = public )

        empty  = DoubleDouble.MakeStatic( [], [] );
        zero   = DoubleDouble.MakeStatic( 0, 0 );
        one    = DoubleDouble.MakeStatic( 1, 0 );
        tiny   = DoubleDouble.MakeStatic( 4.93038065763132e-32, 0 );
        pi     = DoubleDouble.MakeStatic( 3.141592653589793116e+00, 1.224646799147353207e-16 );
        piT2   = DoubleDouble.MakeStatic( 6.283185307179586232e+00,  2.449293598294706414e-16 );
        piD2   = DoubleDouble.MakeStatic( 1.570796326794896558e+00,  6.123233995736766036e-17 );
        piD16  = DoubleDouble.MakeStatic( 1.963495408493620697e-01,  7.654042494670957545e-18 );
        log_2  = DoubleDouble.MakeStatic( 6.931471805599452862e-01,  2.319046813846299558e-17 );
        log_10 = DoubleDouble.MakeStatic( 2.302585092994045901e+00, -2.170756223382249351e-16 );

        NInverseFactorial = 15;
        ExpRescale        = 9;
        LogSteps          = 2;

        InverseFactorial = [
            DoubleDouble.MakeStatic( 1.66666666666666657e-01,   9.25185853854297066e-18 );
            DoubleDouble.MakeStatic( 4.16666666666666644e-02,   2.31296463463574266e-18 );
            DoubleDouble.MakeStatic( 8.33333333333333322e-03,   1.15648231731787138e-19 );
            DoubleDouble.MakeStatic( 1.38888888888888894e-03,  -5.30054395437357706e-20 );
            DoubleDouble.MakeStatic( 1.98412698412698413e-04,   1.72095582934207053e-22 );
            DoubleDouble.MakeStatic( 2.48015873015873016e-05,   2.15119478667758816e-23 );
            DoubleDouble.MakeStatic( 2.75573192239858925e-06,  -1.85839327404647208e-22 );
            DoubleDouble.MakeStatic( 2.75573192239858883e-07,   2.37677146222502973e-23 );
            DoubleDouble.MakeStatic( 2.50521083854417202e-08,  -1.44881407093591197e-24 );
            DoubleDouble.MakeStatic( 2.08767569878681002e-09,  -1.20734505911325997e-25 );
            DoubleDouble.MakeStatic( 1.60590438368216133e-10,   1.25852945887520981e-26 );
            DoubleDouble.MakeStatic( 1.14707455977297245e-11,   2.06555127528307454e-28 );
            DoubleDouble.MakeStatic( 7.64716373181981641e-13,   7.03872877733453001e-30 );
            DoubleDouble.MakeStatic( 4.77947733238738525e-14,   4.39920548583408126e-31 );
            DoubleDouble.MakeStatic( 2.81145725434552060e-15,   1.65088427308614326e-31 );
            ];

        SinTable = [
            DoubleDouble.MakeStatic( 0,  0 );
            DoubleDouble.MakeStatic( 1.950903220161282758e-01,  -7.991079068461731263e-18 );
            DoubleDouble.MakeStatic( 3.826834323650897818e-01,  -1.005077269646158761e-17 );
            DoubleDouble.MakeStatic( 5.555702330196021776e-01,   4.709410940561676821e-17 );
            DoubleDouble.MakeStatic( 7.071067811865475727e-01,  -4.833646656726456726e-17 );
            ];

        CosTable = [
            DoubleDouble.MakeStatic( 1,  0 );
            DoubleDouble.MakeStatic( 9.807852804032304306e-01,   1.854693999782500573e-17 );
            DoubleDouble.MakeStatic( 9.238795325112867385e-01,   1.764504708433667706e-17 );
            DoubleDouble.MakeStatic( 8.314696123025452357e-01,   1.407385698472802389e-18 );
            DoubleDouble.MakeStatic( 7.071067811865475727e-01,  -4.833646656726456726e-17 );
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
            elseif isa( in, 'ExtDouble' )
                [ v.v1, v.v2 ] = ToSumOfDoubles( in );
            else
                v.v1 = double( in );
                v.v2 = zeros( size( in ) );
            end
        end

        function v = Promote( ~, v )
            v = DoubleDouble( v );
        end

        function n = PromotionOrder( ~ )
            n = 1;
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = DoubleDouble.MakeStatic( a1, a2 );
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

        function v = NaN( varargin )
            v = DoubleDouble.MakeStatic( NaN( varargin{:}, 'double' ), NaN( varargin{:}, 'double' ) );
        end

        function v = Inf( varargin )
            v = DoubleDouble.MakeStatic( Inf( varargin{:}, 'double' ), Inf( varargin{:}, 'double' ) );
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

    methods ( Static, Access = { ?BaseDoubleDouble, ?BaseExtDoubleProperties, ?ExtDouble } )

        function v = MakeStatic( a1, a2 )
            v = DoubleDouble;
            v.v1 = double( a1 );
            v.v2 = double( a2 );
        end

    end

end
