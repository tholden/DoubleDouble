classdef ( Abstract ) BaseDoubleDouble

    properties ( SetAccess = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties }, GetAccess = public )

        v1
        v2

    end

    methods ( Abstract )

        v = Promote( a, v )
        n = PromotionOrder( v )

    end

    methods ( Abstract, Access = protected )

        v = Make( z, a1, a2 )

    end

    methods ( Access = protected ) % Not sealed, overridden in QuadDouble.

        function v = Plus( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2 ] = DDPlusDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DDPlusUnderlying( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2 ] = DDPlusUnderlying( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingPlusUnderlying( a, b );
                end
            end
            v = c.Make( x1, x2 );
        end

        function v = Times( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2 ] = DDTimesDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DDTimesUnderlying( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2 ] = DDTimesUnderlying( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingTimesUnderlying( a, b );
                end
            end
            v = c.Make( x1, x2 );
        end

        function v = RDivide( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2 ] = DDDividedByDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2 ] = DDDividedByUnderlying( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2 ] = UnderlyingDividedByDD( a, b.v1, b.v2 );
                else
                    [ x1, x2 ] = UnderlyingDividedByUnderlying( a, b );
                end
            end
            v = c.Make( x1, x2 );
        end

        function v = Normalize( v )
            if all( v.v2 == 0, 'all' )
                v.v1 = Normalize( v.v1 );
                if isempty( v.v1 )
                    v.v2 = [];
                else
                    v.v2 = 0;
                end
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
            end
        end

    end

    methods ( Static, Access = protected )

        function [ a, b, c, Na, Nb, Nc ] = JointPromotion( a, b )
            Na = PromotionOrder( a );
            Nb = PromotionOrder( b );
            if ( Na ~= Nb ) && ( abs( Na - Nb ) < 0.5 )
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
                Na = Na - 1;
                a = a.v1;
            end
            if ( Nb > 0 ) && all( b.v2 == 0, 'all' )
                Nb = Nb - 1;
                b = b.v1;
            end
            if Na < Nc
                a = Promote( c.v1, a );
                Na = max( Na, Nc - 1 );
            end
            if Nb < Nc
                b = Promote( c.v1, b );
                Nb = max( Nb, Nc - 1 );
            end
        end

    end

end
