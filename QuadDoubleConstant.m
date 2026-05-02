classdef QuadDoubleConstant < BaseDoubleDouble & ExtDouble

    properties ( Constant )
        empty = QuadDoubleConstant;
        zero = QuadDoubleConstant;
        one = QuadDoubleConstant;
        tiny = QuadDoubleConstant;
        pi = QuadDoubleConstant;
        piT2 = QuadDoubleConstant;
        piD2 = QuadDoubleConstant;
        piD16 = QuadDoubleConstant;
        log_2 = QuadDoubleConstant;
        log_10 = QuadDoubleConstant;
        NInverseFactorial = QuadDoubleConstant;
        ExpRescale = QuadDoubleConstant;
        LogSteps = QuadDoubleConstant;
        InverseFactorial = QuadDoubleConstant;
        SinTable = QuadDoubleConstant;
        CosTable = QuadDoubleConstant;
    end

    methods

        function v = QuadDoubleConstant( ~, ~ )
        end

        function v = Promote( ~, v )
            v = QuadDouble( v );
        end

        function n = PromotionOrder( ~ )
            n = 2;
        end

    end

    methods ( Access = protected )

        function v = Make( ~, a1, a2 )
            v = QuadDoubleConstant.MakeStatic( a1, a2 );
        end

    end

    methods ( Static, Access = { ?BaseDoubleDouble, ?BaseExtDoubleProperties } )

        function v = MakeStatic( a1, a2 )
            v = QuadDoubleConstant;
            v.v1 = a1;
            v.v2 = a2;
        end

    end

end
