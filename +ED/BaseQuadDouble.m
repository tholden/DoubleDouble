classdef ( Abstract ) BaseQuadDouble < ED.BaseDoubleDouble

    methods

        function v = Normalize( v )
            if any(v.v1.v1 == 0.3404149897495286e-64, 'all')
                disp('BaseQuadDouble.Normalize entered!');
                disp('v.v2 before check is:'); disp(v.v2);
                disp('all(v.v2 == 0) is:'); disp(all(v.v2 == 0, 'all'));
            end
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

    methods ( Hidden )

        function [ s1, s2 ] = UnderlyingPlusUnderlying( a, b )
            [ a, b, c, Na, Nb, Nc, FloorNc ] = ED.BaseQuadDouble.JointPromotion( a, b );
            if Na == Nc
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDPlusQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );

                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDPlusDD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDPlusUnderlying( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b );
                end
            elseif Na >= FloorNc - 1
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDPlusDD( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a.v1, a.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDPlusDD( a.v1, a.v2, b.v1, b.v2 );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                else
                    [ x1, x2, x3, x4 ] = DDPlusUnderlying( a.v1, a.v2, b );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDPlusUnderlying( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDPlusUnderlying( b.v1, b.v2, a );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                else
                    [ x1, x2 ] = UnderlyingPlusUnderlying( a, b );
                    x3 = zeros( size( x1 ), 'like', x1 ); x4 = zeros( size( x1 ), 'like', x1 );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                end
            end
            s1 = c.Make( c.v1.Make( x1, x2 ), c.v1.Make( x3, x4 ) );
            s2 = c.Make( c.v1.Make( x5, x6 ), c.v1.Make( x7, x8 ) );
        end

        function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
            [ a, b, c, Na, Nb, Nc, FloorNc ] = ED.BaseQuadDouble.JointPromotion( a, b );
            if Na == Nc
                z = zeros( size( a.v1.v1 ), 'like', a.v1.v1 );
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2 );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b.v1, b.v2, z, z );
                else
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDTimesQD( a.v1.v1, a.v1.v2, a.v2.v1, a.v2.v2, b, z, z, z );
                end
            elseif Na >= FloorNc - 1
                z = zeros( size( b.v1.v1 ), 'like', b.v1.v1 );
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDTimesQD( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a.v1, a.v2, z, z );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDTimesDDAsQD( a.v1, a.v2, b.v1, b.v2 );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                else
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( a.v1, a.v2, b );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                end
            else
                z = zeros( size( b.v1.v1 ), 'like', b.v1.v1 );
                if Nb == Nc
                    [ x1, x2, x3, x4, x5, x6, x7, x8 ] = QDTimesQD( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a, z, z, z );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDTimesUnderlyingAsQD( b.v1, b.v2, a );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                else
                    [ x1, x2 ] = UnderlyingTimesUnderlying( a, b );
                    x3 = zeros( size( x1 ), 'like', x1 ); x4 = zeros( size( x1 ), 'like', x1 );
                    x5 = zeros( size( x1 ), 'like', x1 ); x6 = zeros( size( x1 ), 'like', x1 );
                    x7 = zeros( size( x1 ), 'like', x1 ); x8 = zeros( size( x1 ), 'like', x1 );
                end
            end
            p1 = c.Make( c.v1.Make( x1, x2 ), c.v1.Make( x3, x4 ) );
            p2 = c.Make( c.v1.Make( x5, x6 ), c.v1.Make( x7, x8 ) );
        end

    end

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
                    [ x1, x2, x3, x4 ] = DDPlusDD( a.v1, a.v2, b.v1, b.v2 );
                else
                    [ x1, x2, x3, x4 ] = DDPlusUnderlying( a.v1, a.v2, b );
                end
            else
                if Nb == Nc
                    [ x1, x2, x3, x4 ] = QDPlusUnderlying( b.v1.v1, b.v1.v2, b.v2.v1, b.v2.v2, a );
                elseif Nb >= FloorNc - 1
                    [ x1, x2, x3, x4 ] = DDPlusUnderlying( b.v1, b.v2, a );
                else
                    [ x1, x2 ] = UnderlyingPlusUnderlying( a, b );
                    x3 = zeros( size( x1 ), 'like', x1 ); x4 = zeros( size( x1 ), 'like', x1 );
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
                    [ x1, x2 ] = UnderlyingTimesUnderlying( a, b );
                    x3 = zeros( size( x1 ), 'like', x1 ); x4 = zeros( size( x1 ), 'like', x1 );
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
