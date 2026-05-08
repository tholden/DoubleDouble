classdef QuadDoubleConstant < ED.BaseDoubleDouble & ED.ExtDouble

    properties ( Constant )

        zero  = ED.QuadDoubleConstant;
        one   = ED.QuadDoubleConstant;
        tiny  = ED.QuadDoubleConstant;
        pi    = ED.QuadDoubleConstant;
        log2  = ED.QuadDoubleConstant;
        log10 = ED.QuadDoubleConstant;

        ExpRescale = ED.QuadDoubleConstant;
        LogSteps   = ED.QuadDoubleConstant;
        SqrtSteps  = ED.QuadDoubleConstant;

        InverseFactorial = ED.QuadDoubleConstant;
        SinTable         = ED.QuadDoubleConstant;
        CosTable         = ED.QuadDoubleConstant;

    end

    methods

        function v = QuadDoubleConstant( ~, ~ )
        end

        function v = Promote( ~, v )
            v = QuadDoubleSlow( v );
        end

        function n = PromotionOrder( ~ )
            n = 2;
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

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = ED.QuadDoubleConstant.MakeStatic( a1, a2 );
        end

    end

    methods ( Static, Access = { ?ED.BaseDoubleDouble, ?ED.BaseExtDoubleProperties } )

        function v = MakeStatic( a1, a2 )
            v = ED.QuadDoubleConstant;
            v.v1 = a1;
            v.v2 = a2;
        end

    end

end
