classdef ( Abstract ) BaseQuadDouble < ED.BaseDoubleDouble

    methods ( Access = protected )

        function v = Plus( a, b )
            [ a, b, c, Na, Nb, Nc, FloorNc ] = ED.BaseQuadDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = QDPlusDD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = QDPlusUnderlying( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b );
                end
            elseif Na >= FloorNc - 1
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDPlusDD( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a.v1, a.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDPlusDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDPlusUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDPlusUnderlying( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDPlusUnderlyingAsQD( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingPlusUnderlying( a, b );
                    x3 = 0;
                    x4 = 0;
                end
            end
            v = c.Make( c.v1.Make( x1, x2 ), c.v1.Make( x3, x4 ) );
        end

        function v = Times( a, b )
            [ a, b, c, Na, Nb, Nc, FloorNc ] = ED.BaseQuadDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = QDTimesDD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = QDTimesUnderlying( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b );
                end
            elseif Na >= FloorNc - 1
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDTimesDD( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a.v1, a.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDTimesDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDTimesUnderlying( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( b.v1, b.v2, a );
                else
                    [ x1, x2, x3, x4 ] = UnderlyingTimesUnderlyingAsQD( a, b );
                end
            end
            v = c.Make( c.v1.Make( x1, x2 ), c.v1.Make( x3, x4 ) );
        end

        function v = RDivide( a, b )
            [ a, b, c, Na, Nb, Nc, FloorNc ] = ED.BaseQuadDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDDividedByQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = QDDividedByDD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = QDDividedByUnderlying( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b );
                end
            elseif Na >= FloorNc - 1
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = DDDividedByQD( a.v1, a.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDDividedByDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDDividedByUnderlyingAsQD( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = UnderlyingDividedByQD( a, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = UnderlyingDividedByDDAsQD( a, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = UnderlyingDividedByUnderlyingAsQD( a, b );
                end
            end
            v = c.Make( c.v1.Make( x1, x2 ), c.v1.Make( x3, x4 ) );
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
            if all( v.v1.v2 == 0, 'all' )
                v.v1.v2 = 0;
            end
            if all( v.v2.v2 == 0, 'all' )
                v.v2.v2 = 0;
            end
            v.v1.v1 = Normalize( v.v1.v1 );
            v.v1.v2 = Normalize( v.v1.v2 );
            v.v2.v1 = Normalize( v.v2.v1 );
            v.v2.v2 = Normalize( v.v2.v2 );
            [ v.v1.v1, v.v1.v2, v.v2.v1, v.v2.v2 ] = QDNormalize( v.v1.v1, v.v1.v2, v.v2.v1, v.v2.v2 );
            if all( v.v1.v2 == 0, 'all' )
                if isempty( v.v1.v1 )
                    v.v1.v2 = [];
                else
                    v.v1.v2 = 0;
                end
            end
            if all( v.v2.v2 == 0, 'all' )
                if isempty( v.v2.v1 )
                    v.v2.v2 = [];
                else
                    v.v2.v2 = 0;
                end
            end
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

        function [ a, b, c, Na, Nb, Nc, FloorNc ] = JointPromotion( a, b )
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
            FloorNc = floor( Nc );
            assert( isa( c, 'ED.BaseQuadDouble' ) );
            if ( Na > 0 ) && all( a.v2 == 0, 'all' )
                a = a.v1;
                Na = PromotionOrder( a );
                if ( Na > 0 ) && all( a.v2 == 0, 'all' )
                    a = a.v1;
                    Na = PromotionOrder( a );
                end
            end
            if ( Nb > 0 ) && all( b.v2 == 0, 'all' )
                b = b.v1;
                Nb = PromotionOrder( b );
                if ( Nb > 0 ) && all( b.v2 == 0, 'all' )
                    b = b.v1;
                    Nb = PromotionOrder( b );
                end
            end
            d = c.v1.v1;
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
