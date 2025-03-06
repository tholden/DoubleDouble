function Results = RunTests()
    % RunTests Run all DoubleDouble tests
    %
    % This function runs the comprehensive test suite for the DoubleDouble class
    % and reports the results.
    %
    % Usage:
    %   Results = RunTests()
    %
    % Returns:
    %   Results - TestResult object containing test outcomes
    
    % Import the necessary testing framework
    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    import matlab.unittest.plugins.TAPPlugin;
    import matlab.unittest.plugins.ToFile;
    
    % Create a test suite from the DoubleDoubleTest class
    Suite = TestSuite.fromClass( ?DoubleDoubleTest );
    
    % Create a test runner with verbose output
    Runner = TestRunner.withTextOutput( 'Verbosity', 3 );
    
    % Add TAP plugin to create a TAP-compatible output file
    TapFile = fullfile( pwd, 'TestResults.tap' );
    Runner.addPlugin( TAPPlugin.producingVersion13( ToFile( TapFile ) ) );
    
    % Run the tests
    Results = Runner.run( Suite );
    
    % Display a summary
    disp( ' ' );
    disp( 'Test Summary:' );
    disp( [ '  ' num2str( sum( [ Results.Passed ] ) ) ' tests passed' ] );
    disp( [ '  ' num2str( sum( [ Results.Failed ] ) ) ' tests failed' ] );
    disp( [ '  ' num2str( sum( [ Results.Incomplete ] ) ) ' tests skipped' ] );
    disp( [ '  Total time: ' num2str( sum( [ Results.Duration ] ) ) ' seconds' ] );
    
    % If any tests failed, display details
    if sum( [ Results.Failed ] ) > 0
        disp( ' ' );
        disp( 'Failed Tests:' );
        FailedIdx = find( [ Results.Failed ] );
        for i = 1 : length( FailedIdx )
            TestIdx = FailedIdx( i );
            disp( [ '  ' Results( TestIdx ).Name ] );
            % Get error message - access structure varies by MATLAB version
            try
                ErrorMsg = Results( TestIdx ).Details.Exception.message;
            catch
                try
                    ErrorMsg = Results( TestIdx ).Exception.message;
                catch
                    ErrorMsg = '';
                end
            end
            if ~isempty( ErrorMsg )
                disp( [ '    ' ErrorMsg ] );
            end
        end
    end
end