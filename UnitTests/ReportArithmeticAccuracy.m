vpa_ = @( x ) vpa( x, 300, BoostPrecision = false );

rng( 'default' );

a = QuadDoubleSlow.rand( 100000, 1 );
b = QuadDoubleSlow.rand( 100000, 1 );

% Indices = unique([ 6956, 9046, 9606, 11026, 11125, 15313, 41086, 45919, 53079, 57708, 66528, 69459, 84018, 87435, 91179, 95916, 59508, 77345, 89690, 55764, 46399, 993, 24754, 4716, 33425, 93494, 94544, 6123 ]);

% a = a(Indices);
% b = b(Indices);

%% double addition test
% Construct two dense double numbers
a_ = double( a );
b_ = double( b );
a_vpa = vpa_( a_ );
b_vpa = vpa_( b_ );
s = vpa_( a_ + b_ );
s_vpa = a_vpa + b_vpa;
err = abs( s - s_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'double         addition       error : \t%g\tIndices: %d\n', max_rel, idx(1) );

% Also test double multiplication
p = vpa_( a_ .* b_ );
p_vpa = a_vpa .* b_vpa;
err = abs( p - p_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'double         multiplication error : \t%g\tIndices: %d\n', max_rel, idx(1) );

%% DoubleDouble addition test
% Construct two dense DoubleDouble numbers
a_ = DoubleDouble( a );
b_ = DoubleDouble( b );
a_vpa = vpa_( a_ );
b_vpa = vpa_( b_ );
s = vpa_( a_ + b_ );
s_vpa = a_vpa + b_vpa;
err = abs( s - s_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'DoubleDouble   addition       error : \t%g\tIndices: %d\n', max_rel, idx(1) );

% Also test DoubleDouble multiplication
p = vpa_( a_ .* b_ );
p_vpa = a_vpa .* b_vpa;
err = abs( p - p_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'DoubleDouble   multiplication error : \t%g\tIndices: %d\n', max_rel, idx(1) );

%% QuadDouble addition test
% Construct two dense QuadDouble numbers
a_ = QuadDouble( a );
b_ = QuadDouble( b );
a_vpa = vpa_( a_ );
b_vpa = vpa_( b_ );
s = vpa_( a_ + b_ );
s_vpa = a_vpa + b_vpa;
err = abs( s - s_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'QuadDouble     addition       error : \t%g\tIndices: %d\n', max_rel, idx(1) );

% Also test QuadDouble multiplication
p = vpa_( a_ .* b_ );
p_vpa = a_vpa .* b_vpa;
err = abs( p - p_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'QuadDouble     multiplication error : \t%g\tIndices: %d\n', max_rel, idx(1) );

%% QuadDoubleSlow addition test
% Construct two dense QuadDoubleSlow numbers
a_ = QuadDoubleSlow( a );
b_ = QuadDoubleSlow( b );
a_vpa = vpa_( a_ );
b_vpa = vpa_( b_ );
s = vpa_( a_ + b_ );
s_vpa = a_vpa + b_vpa;
err = abs( s - s_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'QuadDoubleSlow addition       error : \t%g\tIndices: %d\n', max_rel, idx(1) );

% Also test QuadDoubleSlow multiplication
p = vpa_( a_ .* b_ );
p_vpa = a_vpa .* b_vpa;
err = abs( p - p_vpa );
rel = err ./ max( 1, abs( s_vpa ) );
max_rel = max( rel );
idx = find( rel.' == max_rel, 1 );
fprintf( 'QuadDoubleSlow multiplication error : \t%g\tIndices: %d\n', max_rel, idx(1) );
