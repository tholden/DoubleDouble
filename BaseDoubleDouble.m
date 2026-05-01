classdef ( Abstract ) BaseDoubleDouble

    properties ( SetAccess = protected, GetAccess = public )

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
            [ a, b ] = BaseDoubleDouble.JointPromotion( a, b );
            if PromotionOrder( a ) == PromotionOrder( b )
                [ x1, x2 ] = EDPlusED( a.v1, a.v2, b.v1, b.v2 );
                v = a.Make( x1, x2 );
            elseif PromotionOrder( a ) > PromotionOrder( b )
                [ x1, x2 ] = EDPlusUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                v = a.Make( x1, x2 );
            else % PromotionOrder( b ) > PromotionOrder( a )
                [ x1, x2 ] = EDPlusUnderlying( b.v1, b.v2, Promote( b.v1, a ) );
                v = b.Make( x1, x2 );
            end
        end

        function v = Times( a, b )
            [ a, b ] = BaseDoubleDouble.JointPromotion( a, b );
            if PromotionOrder( a ) == PromotionOrder( b )
                [ x1, x2 ] = EDTimesED( a.v1, a.v2, b.v1, b.v2 );
                v = a.Make( x1, x2 );
            elseif PromotionOrder( a ) > PromotionOrder( b )
                [ x1, x2 ] = EDTimesUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                v = a.Make( x1, x2 );
            else % PromotionOrder( b ) > PromotionOrder( a )
                [ x1, x2 ] = EDTimesUnderlying( b.v1, b.v2, Promote( b.v1, a ) );
                v = b.Make( x1, x2 );
            end
        end

        function v = RDivide( a, b )
            [ a, b ] = BaseDoubleDouble.JointPromotion( a, b );
            if PromotionOrder( a ) == PromotionOrder( b )
                [ x1, x2 ] = EDDividedByED( a.v1, a.v2, b.v1, b.v2 );
                v = a.Make( x1, x2 );
            elseif PromotionOrder( a ) > PromotionOrder( b )
                [ x1, x2 ] = EDDividedByUnderlying( a.v1, a.v2, Promote( a.v1, b ) );
                v = a.Make( x1, x2 );
            else % PromotionOrder( b ) > PromotionOrder( a )
                [ x1, x2 ] = EDDividedByUnderlying( b.v1, b.v2, Promote( b.v1, a ) );
                v = b.Make( x1, x2 );
            end
        end

        function v = Normalize( v )
            v.v1 = Normalize( v.v1 );
            v.v2 = Normalize( v.v2 );
            [ v.v1, v.v2 ] = EDNormalize( v.v1, v.v2 );
        end

    end

    methods ( Static, Access = protected )

        function [ a, b ] = JointPromotion( a, b )
            if ( PromotionOrder( a ) ~= PromotionOrder( b ) ) && ( abs( PromotionOrder( a ) - PromotionOrder( b ) ) < 1 )
                if PromotionOrder( a ) > PromotionOrder( b )
                    b = a.Promote( b );
                else
                    a = b.Promote( a );
                end
            end
        end

    end

end
