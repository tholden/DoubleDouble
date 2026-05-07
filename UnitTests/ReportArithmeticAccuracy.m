vpa_ = @( x ) vpa( x, 300, BoostPrecision = false );

a = OctoDouble.rand(1000,1);
b = OctoDouble.rand(1000,1);

%% QuadDouble addition test
% Construct two dense QuadDouble numbers with 8 doubles of significance
a_oct = QuadDouble( a );
b_oct = QuadDouble( b );
a_oct_vpa = vpa_(a_oct);
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct ./ max( 1, abs(s_oct_vpa) );
fprintf('QuadDouble addition error: %g\n', max( rel_oct ) );

% Also test QuadDouble multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul ./ max( 1, abs(p_oct_vpa) );
fprintf('QuadDouble multiplication error: %g\n', max( rel_oct_mul ) );

%% QuadDoubleSlow addition test
% Construct two dense QuadDoubleSlow numbers with 8 doubles of significance
a_oct = QuadDoubleSlow( a );
b_oct = QuadDoubleSlow( b );
a_oct_vpa = vpa_(a_oct);
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct ./ max( 1, abs(s_oct_vpa) );
fprintf('QuadDoubleSlow addition error: %g\n', max( rel_oct ) );

% Also test QuadDoubleSlow multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul ./ max( 1, abs(p_oct_vpa) );
fprintf('QuadDoubleSlow multiplication error: %g\n', max( rel_oct_mul ) );

%% OctoDouble addition test
% Construct two dense OctoDouble numbers with 8 doubles of significance
a_oct = OctoDouble( a );
b_oct = OctoDouble( b );
a_oct_vpa = vpa_(a_oct);
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct ./ max( 1, abs(s_oct_vpa) );
fprintf('OctoDouble addition error: %g\n', max( rel_oct ) );

% Also test OctoDouble multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul ./ max( 1, abs(p_oct_vpa) );
fprintf('OctoDouble multiplication error: %g\n', max( rel_oct_mul ) );

%% OctoDoubleAlt addition test
% Construct two dense OctoDoubleAlt numbers with 8 doubles of significance
a_oct = OctoDoubleAlt( a );
b_oct = OctoDoubleAlt( b );
a_oct_vpa = vpa_(a_oct);
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct ./ max( 1, abs(s_oct_vpa) );
fprintf('OctoDoubleAlt addition error: %g\n', max( rel_oct ) );

% Also test OctoDoubleAlt multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul ./ max( 1, abs(p_oct_vpa) );
fprintf('OctoDoubleAlt multiplication error: %g\n', max( rel_oct_mul ) );

%% OctoDoubleSlow addition test
% Construct two dense OctoDoubleSlow numbers with 8 doubles of significance
a_oct = OctoDoubleSlow( a );
b_oct = OctoDoubleSlow( b );
a_oct_vpa = vpa_(a_oct);
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct ./ max( 1, abs(s_oct_vpa) );
fprintf('OctoDoubleSlow addition error: %g\n', max( rel_oct ) );

% Also test OctoDoubleSlow multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul ./ max( 1, abs(p_oct_vpa) );
fprintf('OctoDoubleSlow multiplication error: %g\n', max( rel_oct_mul ) );

