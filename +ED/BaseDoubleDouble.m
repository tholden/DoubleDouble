classdef ( Abstract ) BaseDoubleDouble

    properties ( SetAccess = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties }, GetAccess = public )

        v1
        v2
        v3 = 0

    end

    methods ( Abstract )

        v = Promote( a, v )
        n = PromotionOrder( v )

    end

    methods ( Abstract, Access = protected )

        v = Make( z, a1, a2, a3 )

    end

    methods ( Access = protected ) % Not sealed, overridden in QuadDouble.

        function v = Plus( a, b )
            [ a, b, c, Na, Nb, Nc, A3, B3 ] = ED.BaseDoubleDouble.JointPromotionV3( a, b );
            if Na == Nc
                if Nb == Nc
                    [ X1, X2, X3 ] = DDPlusDD( a.v1, a.v2, A3, b.v1, b.v2, B3 );
                else
                    [ X1, X2, X3 ] = DDPlusUnderlying( a.v1, a.v2, A3, b );
                end
            else
                if Nb == Nc
                    [ X1, X2, X3 ] = DDPlusUnderlying( b.v1, b.v2, B3, a );
                else
                    [ X1, X2 ] = UnderlyingPlusUnderlying( a, b );
                    X3 = 0;
                end
            end
            v = c.Make( X1, X2, X3 );
        end

        function v = Times( a, b )
            [ a, b, c, Na, Nb, Nc, A3, B3 ] = ED.BaseDoubleDouble.JointPromotionV3( a, b );
            if Na == Nc
                if Nb == Nc
                    [ X1, X2, X3 ] = DDTimesDD( a.v1, a.v2, A3, b.v1, b.v2, B3 );
                else
                    [ X1, X2, X3 ] = DDTimesUnderlying( a.v1, a.v2, A3, b );
                end
            else
                if Nb == Nc
                    [ X1, X2, X3 ] = DDTimesUnderlying( b.v1, b.v2, B3, a );
                else
                    [ X1, X2 ] = UnderlyingTimesUnderlying( a, b );
                    X3 = 0;
                end
            end
            v = c.Make( X1, X2, X3 );
        end

        function v = RDivide( a, b )
            [ a, b, c, Na, Nb, Nc, A3, B3 ] = ED.BaseDoubleDouble.JointPromotionV3( a, b );
            if Na == Nc
                if Nb == Nc
                    [ X1, X2, X3 ] = DDDividedByDD( a.v1, a.v2, A3, b.v1, b.v2, B3 );
                else
                    [ X1, X2, X3 ] = DDDividedByUnderlying( a.v1, a.v2, A3, b );
                end
            else
                if Nb == Nc
                    [ X1, X2, X3 ] = UnderlyingDividedByDD( a, b.v1, b.v2, B3 );
                else
                    [ X1, X2 ] = UnderlyingDividedByUnderlying( a, b );
                    X3 = 0;
                end
            end
            v = c.Make( X1, X2, X3 );
        end

        function v = Normalize( v )
            if all( v.v2 == 0, 'all' )
                v.v1 = Normalize( v.v1 );
                if isempty( v.v1 )
                    v.v2 = [];
                else
                    v.v2 = 0;
                end
                v.v3 = 0;
                return
            end
            v.v1 = Normalize( v.v1 );
            v.v2 = Normalize( v.v2 );
            [ v.v1, v.v2 ] = DDNormalize( v.v1, v.v2 );
            if all( v.v2 == 0, 'all' )
                if isempty( v.v1 )
                    v.v2 = [];
                else
                    v.v2 = 0;
                end
                v.v3 = 0;
            end
        end

    end

    methods ( Static, Access = protected )

        function [ a, b, c, Na, Nb, Nc ] = JointPromotion( a, b )
            Na = PromotionOrder( a );
            Nb = PromotionOrder( b );
            if ( Na ~= Nb ) && ( floor( Na ) == floor( Nb ) )
                if Na > Nb
                    b = a.Promote( b );
                    Nb = Na;
                else
                    a = b.Promote( a );
                    Na = Nb;
                end
            end
            if Na >= Nb
                c = a;
                Nc = Na;
            else
                c = b;
                Nc = Nb;
            end
            if ( Na > 0 ) && all( a.v2 == 0, 'all' )
                a = a.v1;
                Na = PromotionOrder( a );
            end
            if ( Nb > 0 ) && all( b.v2 == 0, 'all' )
                b = b.v1;
                Nb = PromotionOrder( b );
            end
            d = c.v1;
            Nd = PromotionOrder( d );
            if Na < Nd
                a = d.Promote( a );
                Na = Nd;
            end
            if Nb < Nd
                b = d.Promote( b );
                Nb = Nd;
            end
        end

        function [ a, b, c, Na, Nb, Nc, A3, B3 ] = JointPromotionV3( a, b )
            A3 = ED.BaseDoubleDouble.UnpackV3( a );
            B3 = ED.BaseDoubleDouble.UnpackV3( b );
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
        end

    end

    methods ( Static )

        function V3 = UnpackV3( v )
            if ~isa( v, 'ED.BaseDoubleDouble' ) || all( v.v3 == 0, 'all' )
                V3 = 0;
                return
            end
            V3 = v.v3;
            if isa( V3, 'single' )
                [~, E] = log2( abs( double( v.v2 ) + ( v.v2 == 0 ) ) );
                V3 = pow2( double( V3 ), E );
            end
        end

        function V3Total = UnpackV3Total( v )
            V3Total = 0;
            if ~isa( v, 'ED.BaseDoubleDouble' )
                return
            end
            V3 = ED.BaseDoubleDouble.UnpackV3( v );
            if any( V3 ~= 0, 'all' )
                V3Total = double( V3 );
            end
        end

    end

end
