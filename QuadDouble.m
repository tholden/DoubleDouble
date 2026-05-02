classdef QuadDouble < BaseQuadDouble & ExtDouble & QuadDoublePropertiesMixin

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
            if isa( in, 'QuadDouble' ) || isa( in, 'QuadDoubleSlow' ) || isa( in, 'QuadDoubleConstant' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            elseif isa( in, 'DoubleDouble' )
                v.v1 = in;
                v.v2 = DoubleDouble.zeros( size( in ) );
            elseif isa( in, 'ExtDouble' )
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
            n = 2.8;
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

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDouble.MakeStatic( a1, a2 );
        end

    end

    methods ( Static, Access = { ?BaseDoubleDouble, ?BaseExtDoubleProperties } )

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

end
