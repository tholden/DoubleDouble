classdef QuadDouble < ED.BaseQuadDouble & ED.ExtDouble & ED.QuadDoublePropertiesMixin

    methods

        function v = QuadDouble( in, varargin )
            if nargin == 0
                v.v1 = DoubleDouble.empty;
                v.v2 = [];
                return
            end
            if nargin >= 2
                in = QuadDouble( in );
                for i = 1 : length( varargin )
                    in = Plus( in, varargin{ i } );
                end
            end
            if isa( in, 'QuadDouble' ) || isa( in, 'QuadDoubleSlow' ) || isa( in, 'ED.QuadDoubleConstant' )
                v.v1 = in.v1;
                v.v2 = in.v2;
                if isa( in, 'ED.BaseQuadDouble' )
                    v.v3 = in.v3;
                end
            elseif isa( in, 'DoubleDouble' )
                v.v1 = in;
                if isempty( v.v1 )
                    v.v2 = [];
                else
                    v.v2 = 0;
                end
            elseif isa( in, 'ED.ExtDouble' )
                C = cell( 1, 4 );
                [ C{ : } ] = ToSumOfDoubles( in );
                v.v1 = DoubleDouble.MakeStatic( C{ 1 }, C{ 2 } );
                v.v2 = DoubleDouble.MakeStatic( C{ 3 }, C{ 4 } );
            else
                v.v1 = DoubleDouble( double( in ) );
                if isempty( v.v1 )
                    v.v2 = [];
                else
                    v.v2 = 0;
                end
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

        function v = empty( varargin )
            v = QuadDouble.MakeStatic( DoubleDouble.empty( varargin{:} ), [] );
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

        function v = Make( ~, a1, a2, a3 )
            if nargin < 4
                a3 = 0;
            end
            v = QuadDouble.MakeStatic( a1, a2, a3 );
        end

    end

    methods ( Static, Access = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties } )

        function v = MakeStatic( a1, a2, a3 )
            v = QuadDouble;
            v.v1 = DoubleDouble( a1 );
            if isempty( a1 )
                v.v2 = v.v1;
            elseif all( a2 == 0, 'all' )
                v.v2 = 0;
            else
                v.v2 = DoubleDouble( a2 );
            end
            if nargin >= 3 && ~all( a3 == 0, 'all' )
                v.v3 = a3;
            end
        end

    end

end
