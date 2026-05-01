classdef ( Abstract ) BaseQuadDouble < BaseDoubleDouble

    methods ( Access = protected )

        function v = Plus( a, b )
            % [ a, b ] = BaseDoubleDouble.JointPromotion( a, b );
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

            a = BaseQuadDouble.PromoteStatic( a );
            b = BaseQuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

        function v = Times( a, b )
            a = BaseQuadDouble.PromoteStatic( a );
            b = BaseQuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

        function v = RDivide( a, b )
            a = BaseQuadDouble.PromoteStatic( a );
            b = BaseQuadDouble.PromoteStatic( b );
            [ s0, s1, s2, s3 ] = QDDividedByQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

        function v = Normalize( v )
            [ s0, s1, s2, s3 ] = QDNormalize( v.v1.v1, v.v1.v2, v.v2.v1, v.v2.v2 );
            v = QuadDouble.MakeStatic( DoubleDouble.MakeStatic( s0, s1 ), DoubleDouble.MakeStatic( s2, s3 ) );
        end

    end

    methods ( Static, Access = protected )

        function v = PromoteStatic( a ) % TODO Remove
            v = QuadDouble( a );
        end

    end

end
