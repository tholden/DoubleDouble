classdef QuadDoubleConstant < ED.BaseDoubleDouble & ED.ExtDouble

    properties ( Constant )
        empty = ED.QuadDoubleConstant;
        zero = ED.QuadDoubleConstant;
        one = ED.QuadDoubleConstant;
        tiny = ED.QuadDoubleConstant;
        pi = ED.QuadDoubleConstant;
        piT2 = ED.QuadDoubleConstant;
        piD2 = ED.QuadDoubleConstant;
        piD16 = ED.QuadDoubleConstant;
        log_2 = ED.QuadDoubleConstant;
        log_10 = ED.QuadDoubleConstant;
        NInverseFactorial = ED.QuadDoubleConstant;
        ExpRescale = ED.QuadDoubleConstant;
        LogSteps = ED.QuadDoubleConstant;
        InverseFactorial = ED.QuadDoubleConstant;
        SinTable = ED.QuadDoubleConstant;
        CosTable = ED.QuadDoubleConstant;
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
