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

        function [ a1, a2 ] = Split( a )

            if isreal( a )
                Select = ( a > 1.1079139325602226427e+276 ) | ( a < -1.1079139325602226427e+276 ); % 2^917
                a = Assign( a, TimesPowerOf2( Index( a, Select ), 6.1629758220391547298e-33 ), Select ); % 2^( -107 )
                t1 = 81129638414606681695789005144065.0 * a; % 2^106 + 1
                t2 = t1 - a;
                a1 = t1 - t2;
                a2 = a - a1;
                a1 = Assign( a1, TimesPowerOf2( Index( a1, Select ), 162259276829213363391578010288128.0 ), Select ); % 2^107
                a2 = Assign( a2, TimesPowerOf2( Index( a2, Select ), 162259276829213363391578010288128.0 ), Select ); % 2^107
            else
                [ r1, r2 ] = Split( real( a ) );
                [ i1, i2 ] = Split( imag( a ) );
                a1 = complex( r1, i1 );
                a2 = complex( r2, i2 );
            end

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
