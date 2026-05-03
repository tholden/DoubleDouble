classdef QuadDoubleSlow < ED.BaseDoubleDouble & ED.ExtDouble & ED.QuadDoublePropertiesMixin

    methods

        function v = QuadDoubleSlow( in, varargin )
            if nargin == 0
                v.v1 = DoubleDouble.empty;
                v.v2 = DoubleDouble.empty;
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = QuadDoubleSlow.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDoubleSlow' ) || isa( in, 'QuadDouble' ) || isa( in, 'ED.QuadDoubleConstant' )
                v.v1 = in.v1;
                v.v2 = in.v2;
            elseif isa( in, 'DoubleDouble' )
                v.v1 = in;
                v.v2 = 0;
            elseif isa( in, 'ED.ExtDouble' )
                C = cell( 1, 4 );
                [ C{ : } ] = ToSumOfDoubles( in );
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );
            else
                v.v1 = DoubleDouble( double( in ) );
                v.v2 = 0;
            end
        end

        function v = Promote( ~, v )
            v = QuadDoubleSlow( v );
        end

        function n = PromotionOrder( ~ )
            n = 2.1;
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDoubleSlow.MakeStatic( a1, a2 );
        end

    end

    methods ( Static )

        function v = ones( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( ones( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = zeros( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( zeros( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = eye( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( eye( varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

        function v = NaN( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble.NaN( varargin{:} ), DoubleDouble.NaN( varargin{:} ) );
        end

        function v = Inf( varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble.Inf( varargin{:} ), DoubleDouble.Inf( varargin{:} ) );
        end

        function v = randn( varargin )
            v = ToRandn( QuadDoubleSlow.zeros( varargin{:} ) );
        end

        function v = rand( varargin )
            v = ToRand( QuadDoubleSlow.zeros( varargin{:} ) );
        end

        function v = randi( imax, varargin )
            v = QuadDoubleSlow.MakeStatic( DoubleDouble( randi( imax, varargin{:}, 'double' ) ), DoubleDouble( zeros( varargin{:}, 'double' ) ) );
        end

    end

    methods ( Static, Access = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties } )

        function v = MakeStatic( a1, a2 )
            v = QuadDoubleSlow;
            v.v1 = DoubleDouble( a1 );
            if ~isempty( a2 ) && all( a2 == 0, 'all' )
                v.v2 = 0;
            else
                v.v2 = DoubleDouble( a2 );
            end

        end

    end

end
