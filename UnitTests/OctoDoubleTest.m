classdef OctoDoubleTest < matlab.unittest.TestCase

    % OctoDoubleTest Test suite for OctoDouble class

    properties

        AbsTol = 1e-30;  % Absolute tolerance for comparisons
        RelTol = 1e-15;  % Relative tolerance for comparisons

        % Test data
        SmallValues;
        MediumValues;
        LargeValues;
        ComplexValues;
        MatrixValues;

    end

    methods ( TestMethodSetup )

        function CreateTestData( TestCase )
            % Create test data for use in tests
            TestCase.SmallValues = OctoDouble( [ 1e-10, 2e-10, 3e-10 ] );
            TestCase.MediumValues = OctoDouble( [ 1, 2, 3 ] );
            TestCase.LargeValues = OctoDouble( [ 1e10, 2e10, 3e10 ] );
            TestCase.ComplexValues = OctoDouble( [ 1 + 1i, 2 + 2i, 3 + 3i ] );

            % Create matrix test data
            TestCase.MatrixValues = OctoDouble( [ 1, 2, 3; 4, 5, 6; 7, 8, 9 ] );
        end

    end

    methods ( Test )

        function TestCrossValidationVPA( TestCase )
            % Verify that OctoDouble produces identical results to VPA
            % at full extended precision.
            Data = [ 1.123, -2.456; 3.789, -4.012 ];

            VPAData = vpa( sym(Data, 'f'), 135 );
            QDS = OctoDouble( Data );

            % Test basic arithmetic operations
            ResVPA = VPAData * VPAData + VPAData - ( VPAData / 2 );
            ResQDS = QDS * QDS + QDS - ( QDS / 2 );

            [ v1, v2, v3, v4, v5, v6, v7, v8 ] = ToSumOfDoubles( ResQDS );
            ResQDS_VPA = vpa( v1, 135 ) + vpa( v2, 135 ) + vpa( v3, 135 ) + vpa( v4, 135 ) + vpa( v5, 135 ) + vpa( v6, 135 ) + vpa( v7, 135 ) + vpa( v8, 135 );

            TestCase.verifyEqual( double( ResQDS_VPA - ResVPA ), zeros(size(Data)), 'AbsTol', 1e-120 );

            % Test linear algebra routines (LU decomposition)
            % Not supported for VPA matrices efficiently with same exact outputs (P, L, U may differ).
            % But we can skip it for cross validation with VPA since VPA lu can behave differently.
            % Or just skip LU in TestCrossValidation since it is tested in TestLU.
        end

        function TestCrossValidationQuadDouble( TestCase )
            % Verify that OctoDouble closely matches QuadDouble up to QuadDouble precision.
            Data = [ 1.123, -2.456; 3.789, -4.012 ];

            QD  = QuadDouble( Data );
            OD = OctoDouble( Data );

            % Test basic arithmetic operations
            ResQD  = QD * QD + QD - ( QD / 2 );
            ResOD = OD * OD + OD - ( OD / 2 );

            [ v1, v2, v3, v4, ~, ~, ~, ~ ] = ToSumOfDoubles( ResOD );
            [ q1, q2, q3, q4 ] = ToSumOfDoubles( ResQD );

            TestCase.verifyEqual( v1, q1, 'AbsTol', 1e-60 );
            TestCase.verifyEqual( v2, q2, 'AbsTol', 1e-60 );
            TestCase.verifyEqual( v3, q3, 'AbsTol', 1e-60 );
            TestCase.verifyEqual( v4, q4, 'AbsTol', 1e-60 );

            % Test linear algebra routines (LU decomposition)
            [ L_QD, U_QD, ~ ] = lu( QD );
            [ L_OD, U_OD, ~ ] = lu( OD );

            TestCase.verifyEqual( double( L_OD - L_QD ), zeros(size(Data) ), 'AbsTol', 1e-28 );
            TestCase.verifyEqual( double( U_OD - U_QD ), zeros(size(Data) ), 'AbsTol', 1e-28 );
        end

        % Constructor tests
        function TestConstructorEmpty( TestCase )
            A = OctoDouble();
            [ V1, ~, ~, ~, ~, ~, ~, V8 ] = ToSumOfDoubles( A );
            TestCase.verifyEmpty( V1 );
            TestCase.verifyEmpty( V8 );
        end

        function TestConstructorScalar( TestCase )
            A = OctoDouble( 3.14 );
            [ V1, V2, ~, ~, ~, ~, ~, V8 ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, 3.14 );
            TestCase.verifyTrue( abs( V2 ) < TestCase.AbsTol );
            TestCase.verifyTrue( abs( V8 ) < TestCase.AbsTol );
        end

        function TestConstructorArray( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            [ V1, ~, ~, ~, ~, ~, ~, ~ ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, [ 1, 2, 3 ] );
        end

        function TestConstructorComplex( TestCase )
            A = OctoDouble( 1 + 2i );
            [ V1, ~, ~, ~, ~, ~, ~, ~ ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, 1 + 2i );
        end

        function TestConstructorFromOctoDouble( TestCase )
            A = OctoDouble( 3.14 );
            B = OctoDouble( A );
            TestCase.verifyEqual( double( A ), double( B ) );
        end

        % Basic properties tests
        function TestIsReal( TestCase )
            TestCase.verifyTrue( isreal( TestCase.MediumValues ) );
            TestCase.verifyFalse( isreal( TestCase.ComplexValues ) );
        end

        function TestIsFinite( TestCase )
            A = OctoDouble( [ 1, Inf, NaN ] );
            Expected = [ true, false, false ];
            TestCase.verifyEqual( isfinite( A ), Expected );
        end

        function TestIsInf( TestCase )
            A = OctoDouble( [ 1, Inf, NaN ] );
            Expected = [ false, true, false ];
            TestCase.verifyEqual( isinf( A ), Expected );
        end

        function TestIsNaN( TestCase )
            A = OctoDouble( [ 1, Inf, NaN ] );
            Expected = [ false, false, true ];
            TestCase.verifyEqual( isnan( A ), Expected );
        end

        function TestSize( TestCase )
            TestCase.verifyEqual( size( TestCase.MatrixValues ), [ 3, 3 ] );
            TestCase.verifyEqual( size( TestCase.MediumValues ), [ 1, 3 ] );
        end

        function TestLength( TestCase )
            TestCase.verifyEqual( length( TestCase.MatrixValues ), 3 );
            TestCase.verifyEqual( length( TestCase.MediumValues ), 3 );
        end

        function TestNumel( TestCase )
            TestCase.verifyEqual( numel( TestCase.MatrixValues ), 9 );
            TestCase.verifyEqual( numel( TestCase.MediumValues ), 3 );
        end

        % Conversion tests
        function TestToDouble( TestCase )
            A = OctoDouble( 3.14 );
            TestCase.verifyEqual( double( A ), 3.14 );
        end

        function TestRealImag( TestCase )
            A = TestCase.ComplexValues;
            Re = real( A );
            Im = imag( A );
            TestCase.verifyEqual( double( Re ), [ 1, 2, 3 ] );
            TestCase.verifyEqual( double( Im ), [ 1, 2, 3 ] );
        end

        function TestConj( TestCase )
            A = TestCase.ComplexValues;
            C = conj( A );
            TestCase.verifyEqual( double( real( C ) ), [ 1, 2, 3 ] );
            TestCase.verifyEqual( double( imag( C ) ), [ -1, -2, -3 ] );
        end

        function TestAngle( TestCase )
            A = OctoDouble( 1 + 1i );
            Ang = angle( A );
            TestCase.verifyEqual( double( Ang ), pi / 4, 'RelTol', TestCase.RelTol );
        end

        % Array manipulation tests
        function TestReshape( TestCase )
            A = TestCase.MatrixValues;
            B = reshape( A, 1, 9 );
            TestCase.verifyEqual( size( B ), [ 1, 9 ] );
            TestCase.verifyEqual( double( B ), [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ] );
        end

        function TestRepmat( TestCase )
            A = OctoDouble( [ 1, 2 ] );
            B = repmat( A, 2, 1 );
            TestCase.verifyEqual( size( B ), [ 2, 2 ] );
            TestCase.verifyEqual( double( B ), [ 1, 2; 1, 2 ] );
        end

        function TestDiag( TestCase )
            A = TestCase.MatrixValues;
            D = diag( A );
            TestCase.verifyEqual( double( D ), [ 1; 5; 9 ] );

            V = OctoDouble( [ 1, 2, 3 ] );
            M = diag( V );
            Expected = eye( 3 ) .* [ 1, 2, 3 ];
            TestCase.verifyEqual( double( M ), Expected );
        end

        function TestTrilTriu( TestCase )
            A = TestCase.MatrixValues;
            L = tril( A );
            U = triu( A );

            ExpectedL = [ 1, 0, 0; 4, 5, 0; 7, 8, 9 ];
            ExpectedU = [ 1, 2, 3; 0, 5, 6; 0, 0, 9 ];

            TestCase.verifyEqual( double( L ), ExpectedL );
            TestCase.verifyEqual( double( U ), ExpectedU );
        end

        % Arithmetic operator tests
        function TestPlus( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            B = OctoDouble( [ 4, 5, 6 ] );
            C = A + B;
            TestCase.verifyEqual( double( C ), [ 5, 7, 9 ] );

            % Test with scalar double
            D = A + 10;
            TestCase.verifyEqual( double( D ), [ 11, 12, 13 ] );
        end

        function TestMinus( TestCase )
            A = OctoDouble( [ 5, 7, 9 ] );
            B = OctoDouble( [ 1, 2, 3 ] );
            C = A - B;
            TestCase.verifyEqual( double( C ), [ 4, 5, 6 ] );

            % Test with scalar double
            D = A - 5;
            TestCase.verifyEqual( double( D ), [ 0, 2, 4 ] );
        end

        function TestUminus( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            C = -A;
            TestCase.verifyEqual( double( C ), [ -1, -2, -3 ] );
        end

        function TestTimes( TestCase )
            A = OctoDouble( [ 2, 3, 4 ] );
            B = OctoDouble( [ 5, 6, 7 ] );
            C = A .* B;
            TestCase.verifyEqual( double( C ), [ 10, 18, 28 ] );

            % Test with scalar double
            D = A .* 2;
            TestCase.verifyEqual( double( D ), [ 4, 6, 8 ] );
        end

        function TestMtimes( TestCase )
            A = OctoDouble( [ 1, 2; 3, 4 ] );
            B = OctoDouble( [ 5, 6; 7, 8 ] );
            C = A * B;
            TestCase.verifyEqual( double( C ), [ 19, 22; 43, 50 ] );

            % Test with scalar
            D = A * 2;
            TestCase.verifyEqual( double( D ), [ 2, 4; 6, 8 ] );
        end

        function TestRdivide( TestCase )
            A = OctoDouble( [ 10, 20, 30 ] );
            B = OctoDouble( [ 2, 4, 5 ] );
            C = A ./ B;
            TestCase.verifyEqual( double( C ), [ 5, 5, 6 ] );

            % Test with scalar double
            D = A ./ 2;
            TestCase.verifyEqual( double( D ), [ 5, 10, 15 ] );
        end

        function TestLdivide( TestCase )
            A = OctoDouble( 2 );
            B = OctoDouble( [ 10, 20, 30 ] );
            C = A .\ B;
            TestCase.verifyEqual( double( C ), [ 5, 10, 15 ] );
        end

        function TestMldivide( TestCase )
            A = OctoDouble( [ 3, 1; 1, 2 ] );
            B = OctoDouble( [ 9; 8 ] );
            X = A \ B;
            Expected = [ 2; 3 ];
            TestCase.verifyEqual( double( X ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestPower( TestCase )
            A = OctoDouble( [ 2, 3, 4 ] );
            B = A .^ 2;
            Expected = [ 4, 9, 16 ];
            TestCase.verifyEqual( double( B ), Expected );

            % Test with negative powers
            C = A .^ ( -1 );
            Expected = [ 0.5, 1 / 3, 0.25 ];
            TestCase.verifyEqual( double( C ), Expected, 'RelTol', TestCase.RelTol );

            % Test with fractional powers
            D = A .^ 0.5;
            Expected = sqrt( [ 2, 3, 4 ] );
            TestCase.verifyEqual( double( D ), Expected, 'RelTol', TestCase.RelTol );
        end

        % Comparison operator tests
        function TestComparisons( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            B = OctoDouble( [ 3, 2, 1 ] );

            TestCase.verifyEqual( A < B, [ true, false, false ] );
            TestCase.verifyEqual( A > B, [ false, false, true ] );
            TestCase.verifyEqual( A <= B, [ true, true, false ] );
            TestCase.verifyEqual( A >= B, [ false, true, true ] );
            TestCase.verifyEqual( A == B, [ false, true, false ] );
            TestCase.verifyEqual( A ~= B, [ true, false, true ] );
        end

        % Array construction tests
        function TestColon( TestCase )
            A = OctoDouble( 1 ) : OctoDouble( 5 );
            TestCase.verifyEqual( double( A ), 1 : 5 );

            B = OctoDouble( 1 ) : OctoDouble( 0.5 ) : OctoDouble( 3 );
            TestCase.verifyEqual( double( B ), 1 : 0.5 : 3 );
        end

        function TestConcatenation( TestCase )
            A = OctoDouble( [ 1, 2 ] );
            B = OctoDouble( [ 3, 4 ] );

            C = [ A, B ];
            TestCase.verifyEqual( double( C ), [ 1, 2, 3, 4 ] );

            D = [ A; B ];
            TestCase.verifyEqual( double( D ), [ 1, 2; 3, 4 ] );
        end

        % High-level math operation tests
        function TestSum( TestCase )
            A = TestCase.MatrixValues;
            S1 = sum( A );
            S2 = sum( A, 2 );

            TestCase.verifyEqual( double( S1 ), [ 12, 15, 18 ] );
            TestCase.verifyEqual( double( S2 ), [ 6; 15; 24 ] );
        end

        function TestProd( TestCase )
            A = OctoDouble( [ 2, 3, 4 ] );
            P = prod( A );
            TestCase.verifyEqual( double( P ), 24 );
        end

        function TestCumsum( TestCase )
            A = OctoDouble( [ 1, 2, 3, 4 ] );
            CS = cumsum( A );
            TestCase.verifyEqual( double( CS ), [ 1, 3, 6, 10 ], 'AbsTol', TestCase.AbsTol );
        end

        function TestCumprod( TestCase )
            A = OctoDouble( [ 1, 2, 3, 4 ] );
            CP = cumprod( A );
            TestCase.verifyEqual( double( CP ), [ 1, 2, 6, 24 ], 'AbsTol', TestCase.AbsTol );
        end

        function TestDiff( TestCase )
            A = OctoDouble( [ 1, 3, 6, 10 ] );
            D = diff( A );
            TestCase.verifyEqual( double( D ), [ 2, 3, 4 ], 'AbsTol', TestCase.AbsTol );
        end

        function TestDot( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            B = OctoDouble( [ 4, 5, 6 ] );
            D = dot( A, B );
            TestCase.verifyEqual( double( D ), 32 );
        end

        function TestNorm( TestCase )
            A = OctoDouble( [ 3, 4 ] );
            N = norm( A );
            TestCase.verifyEqual( double( N ), 5 );

            % Test with p-norm
            A = OctoDouble( [ 1, 2, 3 ] );
            N1 = norm( A, 1 );
            NInf = norm( A, Inf );
            TestCase.verifyEqual( double( N1 ), 6 );
            TestCase.verifyEqual( double( NInf ), 3 );
        end

        function TestAbs( TestCase )
            A = OctoDouble( [ -1, 2, -3 ] );
            B = abs( A );
            TestCase.verifyEqual( double( B ), [ 1, 2, 3 ] );

            % Test with complex values
            C = OctoDouble( [ 3+4i, 0, 1-1i ] );
            D = abs( C );
            TestCase.verifyEqual( double( D ), [ 5, 0, sqrt(2) ], 'RelTol', TestCase.RelTol );
        end

        function TestSign( TestCase )
            A = OctoDouble( [ -5, 0, 3 ] );
            S = sign( A );
            TestCase.verifyEqual( double( S ), [ -1, 0, 1 ] );
        end

        % Rounding functions tests
        function TestFloor( TestCase )
            A = OctoDouble( [ 1.7, -1.7 ] );
            F = floor( A );
            TestCase.verifyEqual( double( F ), [ 1, -2 ] );
        end

        function TestCeil( TestCase )
            A = OctoDouble( [ 1.3, -1.3 ] );
            C = ceil( A );
            TestCase.verifyEqual( double( C ), [ 2, -1 ] );
        end

        function TestFix( TestCase )
            A = OctoDouble( [ 1.7, -1.7 ] );
            F = fix( A );
            TestCase.verifyEqual( double( F ), [ 1, -1 ] );
        end

        function TestRound( TestCase )
            A = OctoDouble( [ 1.4, 1.5, 2.5, -1.5 ] );
            R = round( A );
            TestCase.verifyEqual( double( R ), [ 1, 2, 3, -2 ] );
        end

        % Exponential and logarithmic functions tests
        function TestSqrt( TestCase )
            A = OctoDouble( [ 4, 9, 16 ] );
            S = sqrt( A );
            TestCase.verifyEqual( double( S ), [ 2, 3, 4 ] );
        end

        function TestRealsqrt( TestCase )
            A = OctoDouble( [ 4, 9, 16 ] );
            S = realsqrt( A );
            TestCase.verifyEqual( double( S ), [ 2, 3, 4 ] );
        end
        function TestExp( TestCase )
            A = OctoDouble( [ 0, 1, log( 10 ) ] );
            E = exp( A );
            Expected = [ 1, exp( 1 ), 10 ];
            TestCase.verifyEqual( double( E ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestExpm1( TestCase )
            A = OctoDouble( [ 0, 1e-10, 1 ] );
            E = expm1( A );
            Expected = [ 0, 1.0000000000500000000e-10, 1.71828182845904523536028747135 ];
            TestCase.verifyEqual( double( E ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestLog( TestCase )
            A = OctoDouble( [ 1, exp( 1 ), 10 ] );
            L = log( A );
            Expected = [ 0, 1, log( 10 ) ];
            TestCase.verifyEqual( double( L ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestLog10( TestCase )
            A = OctoDouble( [ 1, 10, 100 ] );
            L = log10( A );
            TestCase.verifyEqual( double( L ), [ 0, 1, 2 ], 'RelTol', TestCase.RelTol );
        end

        function TestLog2( TestCase )
            A = OctoDouble( [ 1, 2, 4, 8 ] );
            L = log2( A );
            TestCase.verifyEqual( double( L ), [ 0, 1, 2, 3 ], 'RelTol', TestCase.RelTol );
        end

        % Trigonometric functions tests
        function TestSin( TestCase )
            A = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3, OctoDouble.pi / 2 ];
            S = sin( A );
            Expected = [ 0, 0.5, 1 / sqrt( OctoDouble( 2 ) ), sqrt( OctoDouble( 3 ) ) / 2, 1 ];
            TestCase.verifyEqual( double( S ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestCos( TestCase )
            A = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3, OctoDouble.pi / 2 ];
            C = cos( A );
            Expected = [ 1, sqrt( OctoDouble( 3 ) ) / 2, 1 / sqrt( OctoDouble( 2 ) ), 0.5, 0 ];
            TestCase.verifyEqual( double( C ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestTan( TestCase )
            A = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3 ];
            T = tan( A );
            Expected = [ 0, 1 / sqrt( OctoDouble( 3 ) ), 1, sqrt( OctoDouble( 3 ) ) ];
            TestCase.verifyEqual( double( T ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestAsin( TestCase )
            A = [ 0, 0.5, 1 / sqrt( OctoDouble( 2 ) ), sqrt( OctoDouble( 3 ) ) / 2, 1 ];
            As = asin( A );
            Expected = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3, OctoDouble.pi / 2 ];
            TestCase.verifyEqual( double( As ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestAcos( TestCase )
            A = [ 1, sqrt( OctoDouble( 3 ) ) / 2, 1 / sqrt( OctoDouble( 2 ) ), 0.5, 0 ];
            Ac = acos( A );
            Expected = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3, OctoDouble.pi / 2 ];
            TestCase.verifyEqual( double( Ac ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestAtan( TestCase )
            A = [ 0, 1 / sqrt( OctoDouble( 3 ) ), 1, sqrt( OctoDouble( 3 ) ) ];
            At = atan( A );
            Expected = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3 ];
            TestCase.verifyEqual( double( At ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        function TestAtan2( TestCase )
            Y = OctoDouble( [ 0, 1, 1, 1 ] );
            X = [ 1, sqrt( OctoDouble( 3 ) ), 1, 1 / sqrt( OctoDouble( 3 ) ) ];
            A = atan2( Y, X );
            Expected = [ 0, OctoDouble.pi / 6, OctoDouble.pi / 4, OctoDouble.pi / 3 ];
            TestCase.verifyEqual( double( A ), double( Expected ), 'RelTol', TestCase.RelTol );
        end

        % Hyperbolic functions tests
        function TestSinh( TestCase )
            A = OctoDouble( [ 0, 1, 2 ] );
            S = sinh( A );
            Expected = [ 0, sinh( 1 ), sinh( 2 ) ];
            TestCase.verifyEqual( double( S ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestCosh( TestCase )
            A = OctoDouble( [ 0, 1, 2 ] );
            C = cosh( A );
            Expected = [ 1, cosh( 1 ), cosh( 2 ) ];
            TestCase.verifyEqual( double( C ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestTanh( TestCase )
            A = OctoDouble( [ 0, 1, 2 ] );
            T = tanh( A );
            Expected = [ 0, tanh( 1 ), tanh( 2 ) ];
            TestCase.verifyEqual( double( T ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestAsinh( TestCase )
            A = OctoDouble( [ 0, 1, 2 ] );
            As = asinh( A );
            Expected = [ 0, asinh( 1 ), asinh( 2 ) ];
            TestCase.verifyEqual( double( As ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestAcosh( TestCase )
            A = OctoDouble( [ 1, 2, 3 ] );
            Ac = acosh( A );
            Expected = [ 0, acosh( 2 ), acosh( 3 ) ];
            TestCase.verifyEqual( double( Ac ), Expected, 'RelTol', TestCase.RelTol );
        end

        function TestAtanh( TestCase )
            A = OctoDouble( [ 0, 0.5, 0.75 ] );
            At = atanh( A );
            Expected = [ 0, atanh( 0.5 ), atanh( 0.75 ) ];
            TestCase.verifyEqual( double( At ), Expected, 'RelTol', TestCase.RelTol );
        end

        % Remainder functions tests
        function TestMod( TestCase )
            A = OctoDouble( [ 7, 7, -7, -7 ] );
            B = OctoDouble( [ 3, -3, 3, -3 ] );
            M = mod( A, B );
            Expected = mod( double( A ), double( B ) );
            TestCase.verifyEqual( double( M ), Expected );
        end

        function TestRem( TestCase )
            A = OctoDouble( [ 7, 7, -7, -7 ] );
            B = OctoDouble( [ 3, -3, 3, -3 ] );
            R = rem( A, B );
            Expected = [ 1, 1, -1, -1 ];
            TestCase.verifyEqual( double( R ), Expected );
        end

        % Matrix functions tests
        function TestLU( TestCase )
            A = OctoDouble( [ 2, -1, 0; -1, 2, -1; 0, -1, 2 ] );
            [ L, U, P ] = lu( A );

            % Verify that P*A = L*U
            PA = P * A;
            LU = L * U;
            TestCase.verifyEqual( double( PA ), double( LU ), 'RelTol', TestCase.RelTol );
        end

        function TestQR( TestCase )
            A = OctoDouble( [ 12, -51, 4; 6, 167, -68; -4, 24, -41 ] );
            [ Q, R ] = qr( A );

            % Verify that A = Q*R
            QR = Q * R;
            TestCase.verifyEqual( double( QR ), double( A ), 'RelTol', TestCase.RelTol );

            % Verify that Q is orthogonal
            I = Q' * Q;
            EyeVal = eye( size( I ) );
            TestCase.verifyEqual( double( I ), EyeVal, 'AbsTol', TestCase.AbsTol );
        end

        function TestDet( TestCase )
            A = OctoDouble( [ 1, 2; 3, 4 ] );
            D = det( A );
            TestCase.verifyEqual( double( D ), -2, 'AbsTol', TestCase.AbsTol );

            B = OctoDouble( [ 1, 2, 3; 4, 5, 6; 7, 8, 9 ] );
            D = det( B );
            TestCase.verifyEqual( double( D ), 0, 'AbsTol', TestCase.AbsTol );
        end

        function TestInv( TestCase )
            A = OctoDouble( [ 4, 7; 2, 6 ] );
            AInv = inv( A );

            % Verify that A*AInv = I
            I = A * AInv; %#ok<MINV>
            EyeVal = eye( size( I ) );
            TestCase.verifyEqual( double( I ), EyeVal, 'AbsTol', TestCase.AbsTol, 'RelTol', TestCase.RelTol );
        end

        function TestChol( TestCase )
            A = OctoDouble( [ 4, 12, -16; 12, 37, -43; -16, -43, 98 ] );
            R = chol( A );

            % Verify that R'*R = A
            RtR = R' * R;
            TestCase.verifyEqual( double( RtR ), double( A ), 'RelTol', TestCase.RelTol );
        end

        function TestLDL( TestCase )
            A = OctoDouble( [ 4, 12, -16; 12, 37, -43; -16, -43, 98 ] );
            [ L, D ] = ldl( A, 'vector' );

            % Verify that L*diag(D)*L' = A
            LDLt = L * D * L';
            TestCase.verifyEqual( double( LDLt ), double( A ), 'RelTol', TestCase.RelTol );
        end
        % Tests for newly added functions
        function TestUnique( TestCase )
            A = OctoDouble( [ 2, 1, 2, 3, 1, 4 ] );
            [ C, ia, ic ] = unique( A );
            [ C_, ia_, ic_ ] = unique( double( A ) );

            TestCase.verifyEqual( double( C ), C_ );
            TestCase.verifyEqual( double( ia ), ia_ );
            TestCase.verifyEqual( double( ic ), ic_ );
        end

        function TestMean( TestCase )
            A = OctoDouble( [ 1, 2, 3, 4 ] );
            M = mean( A );
            TestCase.verifyEqual( double( M ), 2.5 );

            % Test with matrix
            B = TestCase.MatrixValues;
            M1 = mean( B, 1 );
            M2 = mean( B, 2 );
            TestCase.verifyEqual( double( M1 ), [ 4, 5, 6 ] );
            TestCase.verifyEqual( double( M2 ), [ 2; 5; 8 ] );
        end

        function TestMedian( TestCase )
            A = OctoDouble( [ 1, 3, 5, 7 ] );
            M = median( A );
            TestCase.verifyEqual( double( M ), 4 );

            B = OctoDouble( [ 1, 3, 5 ] );
            M = median( B );
            TestCase.verifyEqual( double( M ), 3 );
        end

        function TestStd( TestCase )
            A = OctoDouble( [ 1, 2, 3, 4, 5 ] );
            S0 = std( A, 0 );  % Normalize by N - 1
            S1 = std( A, 1 );  % Normalize by N

            Expected0 = std( [ 1, 2, 3, 4, 5 ], 0 );
            Expected1 = std( [ 1, 2, 3, 4, 5 ], 1 );

            TestCase.verifyEqual( double( S0 ), Expected0, 'RelTol', TestCase.RelTol );
            TestCase.verifyEqual( double( S1 ), Expected1, 'RelTol', TestCase.RelTol );
        end

        function TestVar( TestCase )
            A = OctoDouble( [ 1, 2, 3, 4, 5 ] );
            V0 = var( A, 0 );  % Normalize by N - 1
            V1 = var( A, 1 );  % Normalize by N

            Expected0 = var( [ 1, 2, 3, 4, 5 ], 0 );
            Expected1 = var( [ 1, 2, 3, 4, 5 ], 1 );

            TestCase.verifyEqual( double( V0 ), Expected0, 'RelTol', TestCase.RelTol );
            TestCase.verifyEqual( double( V1 ), Expected1, 'RelTol', TestCase.RelTol );
        end

        function TestMeshgrid( TestCase )
            X = OctoDouble( [ 1, 2, 3 ] );
            Y = OctoDouble( [ 4, 5 ] );
            [ XGrid, YGrid ] = meshgrid( X, Y );

            ExpectedX = [ 1, 2, 3; 1, 2, 3 ];
            ExpectedY = [ 4, 4, 4; 5, 5, 5 ];

            TestCase.verifyEqual( double( XGrid ), ExpectedX );
            TestCase.verifyEqual( double( YGrid ), ExpectedY );
        end

        function TestLinspace( TestCase )
            A = OctoDouble( 1 );
            B = OctoDouble( 5 );
            Y = linspace( A, B, 5 );

            Expected = [ 1, 2, 3, 4, 5 ];
            TestCase.verifyEqual( double( Y ), Expected );
        end

        % Test the static methods
        function TestStaticOnes( TestCase )
            A = OctoDouble.ones( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyEqual( double( A ), ones( 2, 3 ) );
        end

        function TestStaticZeros( TestCase )
            A = OctoDouble.zeros( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyEqual( double( A ), zeros( 2, 3 ) );
        end

        function TestStaticEye( TestCase )
            A = OctoDouble.eye( 3 );
            TestCase.verifyEqual( size( A ), [ 3, 3 ] );
            TestCase.verifyEqual( double( A ), eye( 3 ) );
        end

        function TestStaticRand( TestCase )
            A = OctoDouble.rand( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyTrue( all( A >= 0 & A <= 1, 'all' ) );
        end

        function TestStaticRandn( TestCase )
            A = OctoDouble.randn( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
        end

        % Regression tests for bugs found in code review
        function TestCummaxCumminWithDim( TestCase )
            % Bug #2: cummax/cummin ignored the Dim argument
            A = OctoDouble( [ 1, 3; 4, 2 ] );

            % cummax along dim 1 (columns)
            CM1 = cummax( A, 1 );
            TestCase.verifyEqual( double( CM1 ), [ 1, 3; 4, 3 ] );

            % cummax along dim 2 (rows)
            CM2 = cummax( A, 2 );
            TestCase.verifyEqual( double( CM2 ), [ 1, 3; 4, 4 ] );

            % cummin along dim 1 (columns)
            Cm1 = cummin( A, 1 );
            TestCase.verifyEqual( double( Cm1 ), [ 1, 3; 1, 2 ] );

            % cummin along dim 2 (rows)
            Cm2 = cummin( A, 2 );
            TestCase.verifyEqual( double( Cm2 ), [ 1, 1; 4, 2 ] );
        end

        function TestColonMixedTypes( TestCase )
            % Bug #3: colon with plain double start and ED.ExtDouble end
            A = OctoDouble( 1 ) : 5;
            TestCase.verifyEqual( double( A ), 1 : 5 );
            TestCase.verifyTrue( isa( A, 'OctoDouble' ) );

            B = OctoDouble( 1 ) : 0.5 : 3;
            TestCase.verifyEqual( double( B ), 1 : 0.5 : 3 );
            TestCase.verifyTrue( isa( B, 'OctoDouble' ) );
        end

        function TestConv( TestCase )
            % Bug #4: conv crashed due to u.Dot call
            U = OctoDouble( [ 1, 2, 3 ] );
            V = OctoDouble( [ 1, 1 ] );
            W = conv( U, V );
            TestCase.verifyEqual( double( W ), conv( [ 1, 2, 3 ], [ 1, 1 ] ) );
        end

        function TestMrdivideMatrix( TestCase )
            % Bug #5/6: mrdivide for matrices (A / B)
            A = OctoDouble( [ 3, 1; 1, 2 ] );
            B = OctoDouble( [ 1, 0; 0, 1 ] );
            X = A / B;
            TestCase.verifyEqual( double( X ), double( A ), 'RelTol', TestCase.RelTol );

            % Non-trivial case: A / B where B is not identity
            C = OctoDouble( [ 2, 1; 1, 3 ] );
            X2 = A / C;
            % Verify: X2 * C = A
            TestCase.verifyEqual( double( X2 * C ), double( A ), 'RelTol', TestCase.RelTol );
        end

    end

end
