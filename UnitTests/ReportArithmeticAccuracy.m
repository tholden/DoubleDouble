vpa_ = @( x ) vpa( x, 300, BoostPrecision = false );

a = OctoDouble.rand(1000,1);
b = OctoDouble.rand(1000,1);

%% QuadDouble addition test
% Construct two dense QuadDouble numbers with 8 doubles of significance
a_ = QuadDouble( a );
b_ = QuadDouble( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
fprintf('QuadDouble addition error: %g\n', max( rel ) );

% Also test QuadDouble multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err_mul = abs(p - p_vpa);
rel_mul = err_mul ./ max( 1, abs(p_vpa) );
fprintf('QuadDouble multiplication error: %g\n', max( rel_mul ) );

%% QuadDoubleSlow addition test
% Construct two dense QuadDoubleSlow numbers with 8 doubles of significance
a_ = QuadDoubleSlow( a );
b_ = QuadDoubleSlow( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
fprintf('QuadDoubleSlow addition error: %g\n', max( rel ) );

% Also test QuadDoubleSlow multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err_mul = abs(p - p_vpa);
rel_mul = err_mul ./ max( 1, abs(p_vpa) );
fprintf('QuadDoubleSlow multiplication error: %g\n', max( rel_mul ) );

%% OctoDouble addition test
% Construct two dense OctoDouble numbers with 8 doubles of significance
a_ = OctoDouble( a );
b_ = OctoDouble( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
fprintf('OctoDouble addition error: %g\n', max( rel ) );

% Also test OctoDouble multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err_mul = abs(p - p_vpa);
rel_mul = err_mul ./ max( 1, abs(p_vpa) );
fprintf('OctoDouble multiplication error: %g\n', max( rel_mul ) );

%% OctoDoubleAlt addition test
% Construct two dense OctoDoubleAlt numbers with 8 doubles of significance
a_ = OctoDoubleAlt( a );
b_ = OctoDoubleAlt( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
fprintf('OctoDoubleAlt addition error: %g\n', max( rel ) );

% Also test OctoDoubleAlt multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err_mul = abs(p - p_vpa);
rel_mul = err_mul ./ max( 1, abs(p_vpa) );
fprintf('OctoDoubleAlt multiplication error: %g\n', max( rel_mul ) );

%% OctoDoubleSlow addition test
% Construct two dense OctoDoubleSlow numbers with 8 doubles of significance
a_ = OctoDoubleSlow( a );
b_ = OctoDoubleSlow( b );
a_vpa = vpa_(a_);
b_vpa = vpa_(b_);
s = vpa_(a_ + b_);
s_vpa = a_vpa + b_vpa;
err = abs(s - s_vpa);
rel = err ./ max( 1, abs(s_vpa) );
fprintf('OctoDoubleSlow addition error: %g\n', max( rel ) );

% Also test OctoDoubleSlow multiplication
p = vpa_(a_ .* b_);
p_vpa = a_vpa .* b_vpa;
err_mul = abs(p - p_vpa);
rel_mul = err_mul ./ max( 1, abs(p_vpa) );
fprintf('OctoDoubleSlow multiplication error: %g\n', max( rel_mul ) );

