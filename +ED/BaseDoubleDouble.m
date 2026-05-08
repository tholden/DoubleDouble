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
            if Nc == 3
                % OctoDouble output: flatten both operands to 8 raw doubles.
                aa = c.Promote( a );
                bb = c.Promote( b );
                C = cell( 1, 8 );
                [ C{:} ] = ToSumOfDoubles( aa );
                a0 = C{1}; a1_ = C{2}; a2 = C{3}; a3 = C{4}; a4 = C{5}; a5 = C{6}; a6 = C{7}; a7 = C{8};
                [ C{:} ] = ToSumOfDoubles( bb );
                b0 = C{1}; b1_ = C{2}; b2 = C{3}; b3 = C{4}; b4 = C{5}; b5 = C{6}; b6 = C{7}; b7 = C{8};
                [ x0, x1, x2, x3, x4, x5, x6, x7 ] = ODTimesOD( a0, a1_, a2, a3, a4, a5, a6, a7, b0, b1_, b2, b3, b4, b5, b6, b7 );
                v1 = c.v1.Make( c.v1.v1.Make( x0, x1 ), c.v1.v1.Make( x2, x3 ) );
                v2 = c.v1.Make( c.v1.v1.Make( x4, x5 ), c.v1.v1.Make( x6, x7 ) );
                v = c.Make( v1, v2 );
                return
            end
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

    end

end
