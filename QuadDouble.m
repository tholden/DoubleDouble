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

        function v = Normalize( v )
            [ s0, s1, s2, s3 ] = QDNormalize( v.v1.v1, v.v1.v2, v.v2.v1, v.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

    end

    methods ( Static, Access = private )

        function v = PromoteStatic( a ) % TODO Remove
            v = QuadDouble( a );
        end

    end

    methods ( Static, Access = ?ExtDouble )

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
