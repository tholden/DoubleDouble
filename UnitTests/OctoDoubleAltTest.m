classdef OctoDoubleAltTest < matlab.unittest.TestCase

    % OctoDoubleAltTest Test suite for OctoDoubleAlt class

    properties

        Tol = 1e-90;  % Absolute tolerance for comparisons  % Relative tolerance for comparisons

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
            digits( 300 );
            digits( 300 );
            TestCase.SmallValues = OctoDoubleAlt( [ 1e-10, 2e-10, 3e-10 ] );
            TestCase.MediumValues = OctoDoubleAlt( [ 1, 2, 3 ] );
            TestCase.LargeValues = OctoDoubleAlt( [ 1e10, 2e10, 3e10 ] );
            TestCase.ComplexValues = OctoDoubleAlt( [ 1 + 1i, 2 + 2i, 3 + 3i ] );

            % Create matrix test data
            TestCase.MatrixValues = OctoDoubleAlt( [ 1, 2, 3; 4, 5, 6; 7, 8, 9 ] );
        end

    end

    methods ( Test )


        % Constructor tests
        function TestConstructorEmpty( TestCase )
            A = OctoDoubleAlt( );
            [ V1, ~, ~, ~, ~, ~, ~, V8 ] = ToSumOfDoubles( A );
            TestCase.verifyEmpty( V1 );
            TestCase.verifyEmpty( V8 );
        end

        function TestConstructorScalar( TestCase )
            A = OctoDoubleAlt( 3.14 );
            [ V1, V2, ~, ~, ~, ~, ~, V8 ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, 3.14 );
            TestCase.verifyTrue( abs( V2 ) < TestCase.Tol );
            TestCase.verifyTrue( abs( V8 ) < TestCase.Tol );
        end

        function TestConstructorArray( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            [ V1, ~, ~, ~, ~, ~, ~, ~ ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, [ 1, 2, 3 ] );
        end

        function TestConstructorComplex( TestCase )
            A = OctoDoubleAlt( 1 + 2i );
            [ V1, ~, ~, ~, ~, ~, ~, ~ ] = ToSumOfDoubles( A );
            TestCase.verifyEqual( V1, 1 + 2i );
        end

        function TestConstructorFromOctoDouble( TestCase )
            A = OctoDoubleAlt( 3.14 );
            B = OctoDoubleAlt( A );
            TestCase.verifyEqual( double( A ), double( B ) );
        end

        % Basic properties tests
        function TestIsReal( TestCase )
            TestCase.verifyTrue( isreal( TestCase.MediumValues ) );
            TestCase.verifyFalse( isreal( TestCase.ComplexValues ) );
        end

        function TestIsFinite( TestCase )
            A = OctoDoubleAlt( [ 1, Inf, NaN ] );
            Expected = [ true, false, false ];
            TestCase.verifyEqual( isfinite( A ), Expected );
        end

        function TestIsInf( TestCase )
            A = OctoDoubleAlt( [ 1, Inf, NaN ] );
            Expected = [ false, true, false ];
            TestCase.verifyEqual( isinf( A ), Expected );
        end

        function TestIsNaN( TestCase )
            A = OctoDoubleAlt( [ 1, Inf, NaN ] );
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
            A = OctoDoubleAlt( 3.14 );
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
            A = OctoDoubleAlt( 1 + 1i );
            Ang = angle( A );
            expected = vpa( pi ) / 4;
            abs_err = abs( vpa( Ang ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        % Array manipulation tests
        function TestReshape( TestCase )
            A = TestCase.MatrixValues;
            B = reshape( A, 1, 9 );
            TestCase.verifyEqual( size( B ), [ 1, 9 ] );
            TestCase.verifyEqual( double( B ), [ 1, 4, 7, 2, 5, 8, 3, 6, 9 ] );
        end

        function TestRepmat( TestCase )
            A = OctoDoubleAlt( [ 1, 2 ] );
            B = repmat( A, 2, 1 );
            TestCase.verifyEqual( size( B ), [ 2, 2 ] );
            TestCase.verifyEqual( double( B ), [ 1, 2; 1, 2 ] );
        end

        function TestDiag( TestCase )
            A = TestCase.MatrixValues;
            D = diag( A );
            TestCase.verifyEqual( double( D ), [ 1; 5; 9 ] );

            V = OctoDoubleAlt( [ 1, 2, 3 ] );
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
            A = OctoDoubleAlt( [ 1.5, -2.7, 3.1 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            B = OctoDoubleAlt( [ 0.8, 5.5, -1.2 ] ) .* sqrt( OctoDoubleAlt( 2 ) );
            A_vpa = vpa( A );
            B_vpa = vpa( B );

            C = A + B;
            expected = A_vpa + B_vpa;
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            % Test with scalar double
            D = A + 10;
            expectedD = A_vpa + 10;
            abs_errD = abs( vpa( D ) - expectedD );
            errD = abs_errD ./ max( 1, abs( expectedD ) );
            TestCase.verifyLessThanOrEqual( double( max( errD, [], 'all' ) ), TestCase.Tol );
        end


        function TestMinus( TestCase )
            A = OctoDoubleAlt( [ 1.5, -2.7, 3.1 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            B = OctoDoubleAlt( [ 0.8, 5.5, -1.2 ] ) .* sqrt( OctoDoubleAlt( 2 ) );
            A_vpa = vpa( A );
            B_vpa = vpa( B );

            C = A - B;
            expected = A_vpa - B_vpa;
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            % Test with scalar double
            D = A - 5;
            expectedD = A_vpa - 5;
            abs_errD = abs( vpa( D ) - expectedD );
            errD = abs_errD ./ max( 1, abs( expectedD ) );
            TestCase.verifyLessThanOrEqual( double( max( errD, [], 'all' ) ), TestCase.Tol );
        end


        function TestUminus( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            C = -A;
            TestCase.verifyEqual( double( C ), [ -1, -2, -3 ] );
        end

        function TestTimes( TestCase )
            A = OctoDoubleAlt( [ 1.5, -2.7, 3.1 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            B = OctoDoubleAlt( [ 0.8, 5.5, -1.2 ] ) .* sqrt( OctoDoubleAlt( 2 ) );
            A_vpa = vpa( A );
            B_vpa = vpa( B );

            C = A .* B;
            expected = A_vpa .* B_vpa;
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            % Test with scalar double
            D = A .* 2;
            expectedD = A_vpa .* 2;
            abs_errD = abs( vpa( D ) - expectedD );
            errD = abs_errD ./ max( 1, abs( expectedD ) );
            TestCase.verifyLessThanOrEqual( double( max( errD, [], 'all' ) ), TestCase.Tol );
        end


        function TestMtimes( TestCase )
            A = OctoDoubleAlt( [ 1.5, -2.7; 3.1, 0.4 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            B = OctoDoubleAlt( [ 0.8, 5.5; -1.2, 2.3 ] ) .* sqrt( OctoDoubleAlt( 2 ) );
            A_vpa = vpa( A );
            B_vpa = vpa( B );

            C = A * B;
            expected = A_vpa * B_vpa;
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            % Test with scalar
            D = A * 2;
            expectedD = A_vpa * 2;
            abs_errD = abs( vpa( D ) - expectedD );
            errD = abs_errD ./ max( 1, abs( expectedD ) );
            TestCase.verifyLessThanOrEqual( double( max( errD, [], 'all' ) ), TestCase.Tol );
        end


        function TestRdivide( TestCase )
            A = OctoDoubleAlt( [ 1.5, -2.7, 3.1 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            B = OctoDoubleAlt( [ 0.8, 5.5, -1.2 ] ) .* sqrt( OctoDoubleAlt( 2 ) );
            A_vpa = vpa( A );
            B_vpa = vpa( B );

            C = A ./ B;
            expected = A_vpa ./ B_vpa;
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            % Test with scalar double
            D = A ./ 2;
            expectedD = A_vpa ./ 2;
            abs_errD = abs( vpa( D ) - expectedD );
            errD = abs_errD ./ max( 1, abs( expectedD ) );
            TestCase.verifyLessThanOrEqual( double( max( errD, [], 'all' ) ), TestCase.Tol );
        end


        function TestLdivide( TestCase )
            A = OctoDoubleAlt( 2 );
            B = OctoDoubleAlt( [ 10, 20, 30 ] );
            C = A .\ B;
            TestCase.verifyEqual( double( C ), [ 5, 10, 15 ] );
        end

        function TestMldivide( TestCase )
            A = OctoDoubleAlt( [ 3, 1; 1, 2 ] );
            B = OctoDoubleAlt( [ 9; 8 ] );
            X = A \ B;
            Expected = [ 2; 3 ];
            TestCase.verifyEqual( double( X ), Expected, 'RelTol', TestCase.Tol );
        end

        function TestPower( TestCase )
            A = OctoDoubleAlt( [ 2, 3, 4 ] );
            B = A .^ 2;
            Expected = [ 4, 9, 16 ];
            TestCase.verifyEqual( double( B ), Expected );

            % Test with negative powers
            C = A .^ ( -1 );
            Expected = [ 0.5, 1 / 3, 0.25 ];
            TestCase.verifyEqual( double( C ), Expected, 'RelTol', TestCase.Tol );

            % Test with fractional powers
            D = A .^ 0.5;
            Expected = sqrt( [ 2, 3, 4 ] );
            TestCase.verifyEqual( double( D ), Expected, 'RelTol', TestCase.Tol );
        end

        % Comparison operator tests
        function TestComparisons( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            B = OctoDoubleAlt( [ 3, 2, 1 ] );

            TestCase.verifyEqual( A < B, [ true, false, false ] );
            TestCase.verifyEqual( A > B, [ false, false, true ] );
            TestCase.verifyEqual( A <= B, [ true, true, false ] );
            TestCase.verifyEqual( A >= B, [ false, true, true ] );
            TestCase.verifyEqual( A == B, [ false, true, false ] );
            TestCase.verifyEqual( A ~= B, [ true, false, true ] );
        end

        % Array construction tests
        function TestColon( TestCase )
            A = OctoDoubleAlt( 1 ) : OctoDoubleAlt( 5 );
            TestCase.verifyEqual( double( A ), 1 : 5 );

            B = OctoDoubleAlt( 1 ) : OctoDoubleAlt( 0.5 ) : OctoDoubleAlt( 3 );
            TestCase.verifyEqual( double( B ), 1 : 0.5 : 3 );
        end

        function TestConcatenation( TestCase )
            A = OctoDoubleAlt( [ 1, 2 ] );
            B = OctoDoubleAlt( [ 3, 4 ] );

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
            A = OctoDoubleAlt( [ 2, 3, 4 ] );
            P = prod( A );
            TestCase.verifyEqual( double( P ), 24 );
        end

        function TestCumsum( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3, 4 ] );
            CS = cumsum( A );
            TestCase.verifyEqual( double( CS ), [ 1, 3, 6, 10 ], 'AbsTol', TestCase.Tol );
        end

        function TestCumprod( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3, 4 ] );
            CP = cumprod( A );
            TestCase.verifyEqual( double( CP ), [ 1, 2, 6, 24 ], 'AbsTol', TestCase.Tol );
        end

        function TestDiff( TestCase )
            A = OctoDoubleAlt( [ 1, 3, 6, 10 ] );
            D = diff( A );
            TestCase.verifyEqual( double( D ), [ 2, 3, 4 ], 'AbsTol', TestCase.Tol );
        end

        function TestDot( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            B = OctoDoubleAlt( [ 4, 5, 6 ] );
            D = dot( A, B );
            TestCase.verifyEqual( double( D ), 32 );
        end

        function TestNorm( TestCase )
            A = OctoDoubleAlt( [ 3, 4 ] );
            N = norm( A );
            TestCase.verifyEqual( double( N ), 5 );

            % Test with p-norm
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            N1 = norm( A, 1 );
            NInf = norm( A, Inf );
            TestCase.verifyEqual( double( N1 ), 6 );
            TestCase.verifyEqual( double( NInf ), 3 );
        end

        function TestAbs( TestCase )
            A = OctoDoubleAlt( [ -1, 2, -3 ] );
            B = abs( A );
            TestCase.verifyEqual( double( B ), [ 1, 2, 3 ] );

            % Test with complex values
            C = OctoDoubleAlt( [ 3+4i, 0, 1-1i ] );
            D = abs( C );
            TestCase.verifyEqual( double( D ), [ 5, 0, sqrt( 2 ) ], 'RelTol', TestCase.Tol );
        end

        function TestSign( TestCase )
            A = OctoDoubleAlt( [ -5, 0, 3 ] );
            S = sign( A );
            TestCase.verifyEqual( double( S ), [ -1, 0, 1 ] );
        end

        % Rounding functions tests
        function TestFloor( TestCase )
            A = OctoDoubleAlt( [ 1.7, -1.7 ] );
            F = floor( A );
            TestCase.verifyEqual( double( F ), [ 1, -2 ] );
        end

        function TestCeil( TestCase )
            A = OctoDoubleAlt( [ 1.3, -1.3 ] );
            C = ceil( A );
            TestCase.verifyEqual( double( C ), [ 2, -1 ] );
        end

        function TestFix( TestCase )
            A = OctoDoubleAlt( [ 1.7, -1.7 ] );
            F = fix( A );
            TestCase.verifyEqual( double( F ), [ 1, -1 ] );
        end

        function TestRound( TestCase )
            A = OctoDoubleAlt( [ 1.4, 1.5, 2.5, -1.5 ] );
            R = round( A );
            TestCase.verifyEqual( double( R ), [ 1, 2, 3, -2 ] );
        end

        % Exponential and logarithmic functions tests
        function TestSqrt( TestCase )
            A = OctoDoubleAlt( [ 4, 9, 16 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            S = sqrt( A );
            expected = sqrt( vpa( A ) );
            abs_err = abs( vpa( S ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestRealsqrt( TestCase )
            A = OctoDoubleAlt( [ 4, 9, 16 ] ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
            S = realsqrt( A );
            expected = sqrt( vpa( A ) );
            abs_err = abs( vpa( S ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end

        function TestExp( TestCase )
            vals = [0.123, 1.2345, 10.567, 100.23, 0.99] / 10;
            A = OctoDoubleAlt( vals );

            gt_v = vpa( vals, BoostPrecision = false );
            expected = exp( gt_v );

            E = exp( A );
            abs_err = abs( vpa( E ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestExpm1( TestCase )
            A = OctoDoubleAlt( [ 0, 1e-10, 1 ] );
            E = expm1( A );
            Expected = [ 0, 1.0000000000500000000e-10, 1.71828182845904523536028747135 ];
            TestCase.verifyEqual( double( E ), Expected, 'RelTol', TestCase.Tol );
        end

        function TestLog( TestCase )
            vals = [0.123, 1.2345, 10.567, 100.23, 0.99];
            A = OctoDoubleAlt( vals );

            gt_v = vpa( vals, BoostPrecision = false );
            expected = log( gt_v );

            L = log( A );
            abs_err = abs( vpa( L ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestLog10( TestCase )
            A = OctoDoubleAlt( [ 1, 10, 100 ] );
            L = log10( A );
            TestCase.verifyEqual( double( L ), [ 0, 1, 2 ], 'RelTol', TestCase.Tol );
        end

        function TestLog2( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 4, 8 ] );
            L = log2( A );
            TestCase.verifyEqual( double( L ), [ 0, 1, 2, 3 ], 'RelTol', TestCase.Tol );
        end

        % Trigonometric functions tests
        function TestSin( TestCase )
            A = [ 0, OctoDoubleAlt.pi / 6, OctoDoubleAlt.pi / 4, OctoDoubleAlt.pi / 3, OctoDoubleAlt.pi / 2 ];
            S = sin( A );
            expected = sin( vpa( A ) );
            abs_err = abs( vpa( S ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestCos( TestCase )
            A = [ 0, OctoDoubleAlt.pi / 6, OctoDoubleAlt.pi / 4, OctoDoubleAlt.pi / 3, OctoDoubleAlt.pi / 2 ];
            C = cos( A );
            expected = cos( vpa( A ) );
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestTan( TestCase )
            A = [ 0, OctoDoubleAlt.pi / 6, OctoDoubleAlt.pi / 4, OctoDoubleAlt.pi / 3 ];
            T = tan( A );
            expected = tan( vpa( A ) );
            abs_err = abs( vpa( T ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAsin( TestCase )
            A = [ 0, 0.5, 1 / sqrt( OctoDoubleAlt( 2 ) ), sqrt( OctoDoubleAlt( 3 ) ) / 2, 1 ];
            As = asin( A );
            expected = asin( vpa( A ) );
            abs_err = abs( vpa( As ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAcos( TestCase )
            A = [ 1, sqrt( OctoDoubleAlt( 3 ) ) / 2, 1 / sqrt( OctoDoubleAlt( 2 ) ), 0.5, 0 ];
            Ac = acos( A );
            expected = acos( vpa( A ) );
            abs_err = abs( vpa( Ac ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAtan( TestCase )
            A = [ 0, 1 / sqrt( OctoDoubleAlt( 3 ) ), 1, sqrt( OctoDoubleAlt( 3 ) ) ];
            At = atan( A );
            expected = atan( vpa( A ) );
            abs_err = abs( vpa( At ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAtan2( TestCase )
            Y = OctoDoubleAlt( [ 0, 1, 1, 1 ] );
            X = [ 1, sqrt( OctoDoubleAlt( 3 ) ), 1, 1 / sqrt( OctoDoubleAlt( 3 ) ) ];
            A = atan2( Y, X );
            expected = atan2( vpa( Y ), vpa( X ) );
            abs_err = abs( vpa( A ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        % Hyperbolic functions tests
        function TestSinh( TestCase )
            A = OctoDoubleAlt( [ 0, 1, 2 ] );
            S = sinh( A );
            expected = sinh( vpa( A ) );
            abs_err = abs( vpa( S ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestCosh( TestCase )
            A = OctoDoubleAlt( [ 0, 1, 2 ] );
            C = cosh( A );
            expected = cosh( vpa( A ) );
            abs_err = abs( vpa( C ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestTanh( TestCase )
            A = OctoDoubleAlt( [ 0, 1, 2 ] );
            T = tanh( A );
            expected = tanh( vpa( A ) );
            abs_err = abs( vpa( T ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAsinh( TestCase )
            A = OctoDoubleAlt( [ 0, 1, 2 ] );
            As = asinh( A );
            expected = asinh( vpa( A ) );
            abs_err = abs( vpa( As ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAcosh( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3 ] );
            Ac = acosh( A );
            expected = acosh( vpa( A ) );
            abs_err = abs( vpa( Ac ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestAtanh( TestCase )
            A = OctoDoubleAlt( [ 0, 0.5, 0.75 ] );
            At = atanh( A );
            expected = atanh( vpa( A ) );
            abs_err = abs( vpa( At ) - expected );
            err = abs_err ./ max( 1, abs( expected ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        % Remainder functions tests
        function TestMod( TestCase )
            A = OctoDoubleAlt( [ 7, 7, -7, -7 ] );
            B = OctoDoubleAlt( [ 3, -3, 3, -3 ] );
            M = mod( A, B );
            Expected = mod( double( A ), double( B ) );
            TestCase.verifyEqual( double( M ), Expected );
        end

        function TestRem( TestCase )
            A = OctoDoubleAlt( [ 7, 7, -7, -7 ] );
            B = OctoDoubleAlt( [ 3, -3, 3, -3 ] );
            R = rem( A, B );
            Expected = [ 1, 1, -1, -1 ];
            TestCase.verifyEqual( double( R ), Expected );
        end

        % Matrix functions tests
        function TestLU( TestCase )
            A = OctoDoubleAlt( [ 2, -1, 0; -1, 2, -1; 0, -1, 2 ] );
            [ L, U, P ] = lu( A );
            PA = P * A;
            LU = L * U;
            abs_err = abs( vpa( PA ) - vpa( LU ) );
            err = abs_err ./ max( 1, abs( vpa( PA ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestQR( TestCase )
            A = OctoDoubleAlt( [ 12, -51, 4; 6, 167, -68; -4, 24, -41 ] );
            [ Q, R ] = qr( A );
            QR = Q * R;
            abs_err = abs( vpa( QR ) - vpa( A ) );
            err = abs_err ./ max( 1, abs( vpa( A ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
            I = Q' * Q;
            EyeVal = eye( size( I ) );
            abs_err2 = abs( vpa( I ) - vpa( EyeVal ) );
            err2 = abs_err2 ./ max( 1, abs( vpa( EyeVal ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err2, [], 'all' ) ), TestCase.Tol );
        end


        function TestDet( TestCase )
            A = OctoDoubleAlt( [ 1, 2; 3, 4 ] );
            D = det( A );
            expected1 = vpa( -2 );
            abs_err = abs( vpa( D ) - expected1 );
            err = abs_err ./ max( 1, abs( expected1 ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );

            B = OctoDoubleAlt( [ 1, 2, 3; 4, 5, 6; 7, 8, 9 ] );
            D2 = det( B );
            expected2 = vpa( 0 );
            abs_err2 = abs( vpa( D2 ) - expected2 );
            err2 = abs_err2 ./ max( 1, abs( expected2 ) );
            TestCase.verifyLessThanOrEqual( double( max( err2, [], 'all' ) ), TestCase.Tol );
        end


        function TestInv( TestCase )
            A = OctoDoubleAlt( [ 4, 7; 2, 6 ] );
            AInv = inv( A );
            I = A * AInv; %#ok<MINV>
            EyeVal = eye( size( I ) );
            abs_err = abs( vpa( I ) - vpa( EyeVal ) );
            err = abs_err ./ max( 1, abs( vpa( EyeVal ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestChol( TestCase )
            A = OctoDoubleAlt( [ 4, 12, -16; 12, 37, -43; -16, -43, 98 ] );
            R = chol( A );
            RtR = R' * R;
            abs_err = abs( vpa( RtR ) - vpa( A ) );
            err = abs_err ./ max( 1, abs( vpa( A ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end


        function TestLDL( TestCase )
            A = OctoDoubleAlt( [ 4, 12, -16; 12, 37, -43; -16, -43, 98 ] );
            [ L, D ] = ldl( A, 'vector' );
            LDLt = L * D * L';
            abs_err = abs( vpa( LDLt ) - vpa( A ) );
            err = abs_err ./ max( 1, abs( vpa( A ) ) );
            TestCase.verifyLessThanOrEqual( double( max( err, [], 'all' ) ), TestCase.Tol );
        end

        % Tests for newly added functions
        function TestUnique( TestCase )
            A = OctoDoubleAlt( [ 2, 1, 2, 3, 1, 4 ] );
            [ C, ia, ic ] = unique( A );
            [ C_, ia_, ic_ ] = unique( double( A ) );

            TestCase.verifyEqual( double( C ), C_ );
            TestCase.verifyEqual( double( ia ), ia_ );
            TestCase.verifyEqual( double( ic ), ic_ );
        end

        function TestMean( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3, 4 ] );
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
            A = OctoDoubleAlt( [ 1, 3, 5, 7 ] );
            M = median( A );
            TestCase.verifyEqual( double( M ), 4 );

            B = OctoDoubleAlt( [ 1, 3, 5 ] );
            M = median( B );
            TestCase.verifyEqual( double( M ), 3 );
        end

        function TestStd( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3, 4, 5 ] );
            S0 = std( A, 0 );  % Normalize by N - 1
            S1 = std( A, 1 );  % Normalize by N

            Expected0 = std( [ 1, 2, 3, 4, 5 ], 0 );
            Expected1 = std( [ 1, 2, 3, 4, 5 ], 1 );

            TestCase.verifyEqual( double( S0 ), Expected0, 'RelTol', TestCase.Tol );
            TestCase.verifyEqual( double( S1 ), Expected1, 'RelTol', TestCase.Tol );
        end

        function TestVar( TestCase )
            A = OctoDoubleAlt( [ 1, 2, 3, 4, 5 ] );
            V0 = var( A, 0 );  % Normalize by N - 1
            V1 = var( A, 1 );  % Normalize by N

            Expected0 = var( [ 1, 2, 3, 4, 5 ], 0 );
            Expected1 = var( [ 1, 2, 3, 4, 5 ], 1 );

            TestCase.verifyEqual( double( V0 ), Expected0, 'RelTol', TestCase.Tol );
            TestCase.verifyEqual( double( V1 ), Expected1, 'RelTol', TestCase.Tol );
        end

        function TestMeshgrid( TestCase )
            X = OctoDoubleAlt( [ 1, 2, 3 ] );
            Y = OctoDoubleAlt( [ 4, 5 ] );
            [ XGrid, YGrid ] = meshgrid( X, Y );

            ExpectedX = [ 1, 2, 3; 1, 2, 3 ];
            ExpectedY = [ 4, 4, 4; 5, 5, 5 ];

            TestCase.verifyEqual( double( XGrid ), ExpectedX );
            TestCase.verifyEqual( double( YGrid ), ExpectedY );
        end

        function TestLinspace( TestCase )
            A = OctoDoubleAlt( 1 );
            B = OctoDoubleAlt( 5 );
            Y = linspace( A, B, 5 );

            Expected = [ 1, 2, 3, 4, 5 ];
            TestCase.verifyEqual( double( Y ), Expected );
        end

        % Test the static methods
        function TestStaticOnes( TestCase )
            A = OctoDoubleAlt.ones( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyEqual( double( A ), ones( 2, 3 ) );
        end

        function TestStaticZeros( TestCase )
            A = OctoDoubleAlt.zeros( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyEqual( double( A ), zeros( 2, 3 ) );
        end

        function TestStaticEye( TestCase )
            A = OctoDoubleAlt.eye( 3 );
            TestCase.verifyEqual( size( A ), [ 3, 3 ] );
            TestCase.verifyEqual( double( A ), eye( 3 ) );
        end

        function TestStaticRand( TestCase )
            A = OctoDoubleAlt.rand( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
            TestCase.verifyTrue( all( A >= 0 & A <= 1, 'all' ) );
        end

        function TestStaticRandn( TestCase )
            A = OctoDoubleAlt.randn( 2, 3 );
            TestCase.verifyEqual( size( A ), [ 2, 3 ] );
        end

        % Regression tests for bugs found in code review
        function TestCummaxCumminWithDim( TestCase )
            % Bug #2: cummax/cummin ignored the Dim argument
            A = OctoDoubleAlt( [ 1, 3; 4, 2 ] );

            % cummax along dim 1 ( columns )
            CM1 = cummax( A, 1 );
            TestCase.verifyEqual( double( CM1 ), [ 1, 3; 4, 3 ] );

            % cummax along dim 2 ( rows )
            CM2 = cummax( A, 2 );
            TestCase.verifyEqual( double( CM2 ), [ 1, 3; 4, 4 ] );

            % cummin along dim 1 ( columns )
            Cm1 = cummin( A, 1 );
            TestCase.verifyEqual( double( Cm1 ), [ 1, 3; 1, 2 ] );

            % cummin along dim 2 ( rows )
            Cm2 = cummin( A, 2 );
            TestCase.verifyEqual( double( Cm2 ), [ 1, 1; 4, 2 ] );
        end

        function TestColonMixedTypes( TestCase )
            % Bug #3: colon with plain double start and ED.ExtDouble end
            A = OctoDoubleAlt( 1 ) : 5;
            TestCase.verifyEqual( double( A ), 1 : 5 );
            TestCase.verifyTrue( isa( A, 'OctoDoubleAlt' ) );

            B = OctoDoubleAlt( 1 ) : 0.5 : 3;
            TestCase.verifyEqual( double( B ), 1 : 0.5 : 3 );
            TestCase.verifyTrue( isa( B, 'OctoDoubleAlt' ) );
        end

        function TestConv( TestCase )
            % Bug #4: conv crashed due to u.Dot call
            U = OctoDoubleAlt( [ 1, 2, 3 ] );
            V = OctoDoubleAlt( [ 1, 1 ] );
            W = conv( U, V );
            TestCase.verifyEqual( double( W ), conv( [ 1, 2, 3 ], [ 1, 1 ] ) );
        end

        function TestMrdivideMatrix( TestCase )
            % Bug #5/6: mrdivide for matrices ( A / B )
            A = OctoDoubleAlt( [ 3, 1; 1, 2 ] );
            B = OctoDoubleAlt( [ 1, 0; 0, 1 ] );
            X = A / B;
            TestCase.verifyEqual( double( X ), double( A ), 'RelTol', TestCase.Tol );

            % Non-trivial case: A / B where B is not identity
            C = OctoDoubleAlt( [ 2, 1; 1, 3 ] );
            X2 = A / C;
            % Verify: X2 * C = A
            TestCase.verifyEqual( double( X2 * C ), double( A ), 'RelTol', TestCase.Tol );
        end

    end

end
