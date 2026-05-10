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

    methods

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

    methods ( Hidden )

        function [ s1, s2 ] = DDNormalize( a, b )
            [ s1, s2 ] = UnderlyingPlusUnderlying( a, b );
        end

        function [ s1, s2 ] = UnderlyingPlusUnderlying( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDPlusDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDPlusUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDPlusUnderlyingAsQD( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingPlusUnderlying( a, b );
                end
            end
            if exist('x3', 'var')
                s1 = c.Make( x1, x2 );
                s2 = c.Make( x3, x4 );
            else
                s1 = c.Make( x1, x2 );
                s2 = c.zeros( size( s1 ) );
            end
        end

        function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDTimesDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingTimesUnderlying( a, b );
                end
            end
            if exist('x3', 'var')
                p1 = c.Make( x1, x2 );
                p2 = c.Make( x3, x4 );
            else
                p1 = c.Make( x1, x2 );
                p2 = c.zeros( size( p1 ) );
            end
        end

        function [ r1, r2 ] = UnderlyingDividedByUnderlying( a, b )
            [ a, b, c, Na, Nb, Nc ] = ED.BaseDoubleDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDDividedByDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDDividedByUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = UnderlyingDividedByDDAsQD( a, b.v1, b.v2 );
                else
                    [ x1, x2 ] = UnderlyingDividedByUnderlying( a, b );
                end
            end
            if exist('x3', 'var')
                r1 = c.Make( x1, x2 );
                r2 = c.Make( x3, x4 );
            else
                r1 = c.Make( x1, x2 );
                r2 = c.zeros( size( r1 ) );
            end
        end

    end

    methods ( Access = protected )

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

        function v = Minus( a, b )
            v = Plus( a, -b );
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
