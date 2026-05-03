classdef QuadDouble < ED.BaseQuadDouble & ED.ExtDouble & ED.QuadDoublePropertiesMixin

    methods

        function v = QuadDouble( in, varargin )
            if nargin == 0 || isempty( in )
                v.v1 = DoubleDouble.empty;
                v.v2 = DoubleDouble.empty;
                return
            end
            if nargin >= 2
                for i = 1 : length( varargin )
                    in = QuadDouble.Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDouble' ) || isa( in, 'QuadDoubleSlow' ) || isa( in, 'ED.QuadDoubleConstant' )
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
            v = QuadDouble( v );
        end

        function n = PromotionOrder( ~ )
            n = 2.2;
        end

    end

    methods ( Static )

        function v = ones( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.ones( varargin{:} ), 0 );
        end

        function v = zeros( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.zeros( varargin{:} ), 0 );
        end

        function v = eye( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.eye( varargin{:} ), 0 );
        end

        function v = NaN( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.NaN( varargin{:} ), 0 );
        end

        function v = Inf( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.Inf( varargin{:} ), 0 );
        end

        function v = rand( varargin )
            v = ToRand( QuadDouble.zeros( varargin{:} ) );
        end

        function v = randn( varargin )
            v = ToRandn( QuadDouble.zeros( varargin{:} ) );
        end

        function v = randi( imax, varargin )
            v = QuadDouble.MakeStatic( DoubleDouble( randi( imax, varargin{:}, 'double' ) ), 0 );
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDouble.MakeStatic( a1, a2 );
        end

    end

    methods ( Static, Access = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties } )

        function v = MakeStatic( a1, a2 )
            v = QuadDouble;
            v.v1 = DoubleDouble( a1 );
            if isempty( a1 )
                v.v2 = v.v1;
            elseif all( a2 == 0, 'all' )
                v.v2 = 0;
            else
                v.v2 = DoubleDouble( a2 );
            end
        end

    end

end
