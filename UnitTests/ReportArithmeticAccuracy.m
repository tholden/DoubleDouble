vpa_ = @( x ) vpa( x, 300, BoostPrecision = false );

rng( 'default' );

a = OctoDouble.rand(100000,1);
b = OctoDouble.rand(100000,1);

% Indices = [ 6956, 9046, 9606, 11026, 11125, 15313, 41086, 45919, 53079, 57708, 66528, 69459, 84018, 87435, 91179, 95916 ];

% a = a(Indices);
% b = b(Indices);

%% double addition test
% Construct two dense double numbers
a_ = double( a );
b_ = double( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('double         addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test double multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('double         multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% DoubleDouble addition test
% Construct two dense DoubleDouble numbers
a_ = DoubleDouble( a );
b_ = DoubleDouble( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('DoubleDouble   addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test DoubleDouble multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('DoubleDouble   multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% QuadDouble addition test
% Construct two dense QuadDouble numbers
a_ = QuadDouble( a );
b_ = QuadDouble( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('QuadDouble     addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test QuadDouble multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('QuadDouble     multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% QuadDoubleSlow addition test
% Construct two dense QuadDoubleSlow numbers
a_ = QuadDoubleSlow( a );
b_ = QuadDoubleSlow( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('QuadDoubleSlow addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test QuadDoubleSlow multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('QuadDoubleSlow multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% OctoDouble addition test
% Construct two dense OctoDouble numbers
a_ = OctoDouble( a );
b_ = OctoDouble( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDouble     addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test OctoDouble multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDouble     multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% OctoDoubleAlt addition test
% Construct two dense OctoDoubleAlt numbers
a_ = OctoDoubleAlt( a );
b_ = OctoDoubleAlt( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDoubleAlt  addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test OctoDoubleAlt multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDoubleAlt  multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

%% OctoDoubleSlow addition test
% Construct two dense OctoDoubleSlow numbers
a_ = OctoDoubleSlow( a );
b_ = OctoDoubleSlow( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDoubleSlow addition       error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

% Also test OctoDoubleSlow multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err = abs(p - p_vpa);
rel = err ./ max( 1, abs(s_vpa) );
max_rel = max( rel );
fprintf('OctoDoubleSlow multiplication error:\t%g\tIndices: %d\n', max_rel, find( rel.' == max_rel ) );

