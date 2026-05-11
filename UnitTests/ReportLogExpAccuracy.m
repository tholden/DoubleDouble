function TestLogExp()
    addpath( [ pwd filesep '..' ] );

    vals = [ 0.123; 1.2345; 10.567; 100.23; 0.99 ];

    digits( 300 );

    disp( '========================================' );
    disp( '==== Testing Log ====' );
    disp( '========================================' );
    testLog( 'DoubleDouble', vals );
    testLog( 'QuadDoubleSlow', vals );
    testLog( 'QuadDouble', vals );

    disp( '========================================' );
    disp( '==== Testing Exp ====' );
    disp( '========================================' );
    testExp( 'DoubleDouble', vals );
    testExp( 'QuadDoubleSlow', vals );
    testExp( 'QuadDouble', vals );
end

function testLog( className, vals )
    fprintf( ' \ n--- %s ---\n', className);

    % Compute Ground Truth
    gt_v = vpa( vals, BoostPrecision = false );
    gt_log = log( gt_v ); % High precision ground truth

    v = feval( className, vals );

    my_log = vpa( log( v ) );

    % Compute max RELATIVE error across all test values
    abs_err = abs( my_log - gt_log );
    rel_err = abs_err ./abs( gt_log );

    fprintf( 'Abs Error = %g, Rel Error = %g\n', double(max(abs_err)), double(max(rel_err)));
end

function testExp( className, vals )
    fprintf( ' \ n--- %s ---\n', className);

    vals = vals / 10;

    % Compute Ground Truth
    gt_v = vpa( vals, BoostPrecision = false );
    gt_exp = exp( gt_v ); % High precision ground truth

    v = feval( className, vals );

    my_e = vpa( exp( v ) );

    % Compute max RELATIVE error across all test values
    abs_err = abs( my_e - gt_exp );
    rel_err = abs_err ./ abs( gt_exp );

    fprintf( 'Abs Error = %g, Rel Error = %g\n', double(max(abs_err)), double(max(rel_err)));
end
