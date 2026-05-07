vpa_ = @( x ) vpa( x, 300, BoostPrecision = false );

%% QuadDouble addition test
% Construct two dense QuadDouble numbers with 8 doubles of significance
a_oct = QuadDouble( 1.5 ) ./ 7 + exp( QuadDouble( 0.1 ) );
b_oct = QuadDouble( 5.5 ) .* sqrt( QuadDouble( 2 ) );
a_oct_vpa = vpa_(a_oct);
% disp( ' ' );
% disp( 'Check that sum of 8 doubles is exact for QuadDouble:' );
% disp( a_oct_vpa - (  vpa_(1.552722130007543e+03) + vpa_(9.946561756900226e-14) ...
%     + vpa_(2.521714853467526e-30) + (-1.366223618221636e-46)+ vpa_(3.141592653589793e-63) + (-1.234567890123457e-79) ...
%     + vpa_(5.678901234567890e-96) + vpa_(9.876543210987654e-112) ) );
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct / abs(s_oct_vpa);
fprintf('QuadDouble addition absolute error: %g\n', err_oct);
fprintf('QuadDouble addition relative error: %g\n', rel_oct);

% Also test QuadDouble multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul / abs(p_oct_vpa);
fprintf('QuadDouble multiplication absolute error: %g\n', err_oct_mul);
fprintf('QuadDouble multiplication relative error: %g\n', rel_oct_mul);

%% QuadDoubleSlow addition test
% Construct two dense QuadDoubleSlow numbers with 8 doubles of significance
a_oct = QuadDoubleSlow( 1.5 ) ./ 7 + exp( QuadDoubleSlow( 0.1 ) );
b_oct = QuadDoubleSlow( 5.5 ) .* sqrt( QuadDoubleSlow( 2 ) );
a_oct_vpa = vpa_(a_oct);
% disp( ' ' );
% disp( 'Check that sum of 8 doubles is exact for QuadDoubleSlow:' );
% disp( a_oct_vpa - (  vpa_(1.552722130007543e+03) + vpa_(9.946561756900226e-14) ...
%     + vpa_(2.521714853467526e-30) + (-1.366223618221636e-46)+ vpa_(3.141592653589793e-63) + (-1.234567890123457e-79) ...
%     + vpa_(5.678901234567890e-96) + vpa_(9.876543210987654e-112) ) );
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct / abs(s_oct_vpa);
fprintf('QuadDoubleSlow addition absolute error: %g\n', err_oct);
fprintf('QuadDoubleSlow addition relative error: %g\n', rel_oct);

% Also test QuadDoubleSlow multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul / abs(p_oct_vpa);
fprintf('QuadDoubleSlow multiplication absolute error: %g\n', err_oct_mul);
fprintf('QuadDoubleSlow multiplication relative error: %g\n', rel_oct_mul);

%% OctoDouble addition test
% Construct two dense OctoDouble numbers with 8 doubles of significance
a_oct = OctoDouble( 1.5 ) ./ 7 + exp( OctoDouble( 0.1 ) );
b_oct = OctoDouble( 5.5 ) .* sqrt( OctoDouble( 2 ) );
a_oct_vpa = vpa_(a_oct);
% disp( ' ' );
% disp( 'Check that sum of 8 doubles is exact for OctoDouble:' );
% disp( a_oct_vpa - (  vpa_(1.552722130007543e+03) + vpa_(9.946561756900226e-14) ...
%     + vpa_(2.521714853467526e-30) + (-1.366223618221636e-46)+ vpa_(3.141592653589793e-63) + (-1.234567890123457e-79) ...
%     + vpa_(5.678901234567890e-96) + vpa_(9.876543210987654e-112) ) );
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct / abs(s_oct_vpa);
fprintf('OctoDouble addition absolute error: %g\n', err_oct);
fprintf('OctoDouble addition relative error: %g\n', rel_oct);

% Also test OctoDouble multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul / abs(p_oct_vpa);
fprintf('OctoDouble multiplication absolute error: %g\n', err_oct_mul);
fprintf('OctoDouble multiplication relative error: %g\n', rel_oct_mul);

%% OctoDoubleAlt addition test
% Construct two dense OctoDoubleAlt numbers with 8 doubles of significance
a_oct = OctoDoubleAlt( 1.5 ) ./ 7 + exp( OctoDoubleAlt( 0.1 ) );
b_oct = OctoDoubleAlt( 5.5 ) .* sqrt( OctoDoubleAlt( 2 ) );
a_oct_vpa = vpa_(a_oct);
% disp( ' ' );
% disp( 'Check that sum of 8 doubles is exact for OctoDoubleAlt:' );
% disp( a_oct_vpa - (  vpa_(1.552722130007543e+03) + vpa_(9.946561756900226e-14) ...
%     + vpa_(2.521714853467526e-30) + (-1.366223618221636e-46)+ vpa_(3.141592653589793e-63) + (-1.234567890123457e-79) ...
%     + vpa_(5.678901234567890e-96) + vpa_(9.876543210987654e-112) ) );
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct / abs(s_oct_vpa);
fprintf('OctoDoubleAlt addition absolute error: %g\n', err_oct);
fprintf('OctoDoubleAlt addition relative error: %g\n', rel_oct);

% Also test OctoDoubleAlt multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul / abs(p_oct_vpa);
fprintf('OctoDoubleAlt multiplication absolute error: %g\n', err_oct_mul);
fprintf('OctoDoubleAlt multiplication relative error: %g\n', rel_oct_mul);

%% OctoDoubleSlow addition test
% Construct two dense OctoDoubleSlow numbers with 8 doubles of significance
a_oct = OctoDoubleSlow( 1.5 ) ./ 7 + exp( OctoDoubleSlow( 0.1 ) );
b_oct = OctoDoubleSlow( 5.5 ) .* sqrt( OctoDoubleSlow( 2 ) );
a_oct_vpa = vpa_(a_oct);
% disp( ' ' );
% disp( 'Check that sum of 8 doubles is exact for OctoDoubleSlow:' );
% disp( a_oct_vpa - (  vpa_(1.552722130007543e+03) + vpa_(9.946561756900226e-14) ...
%     + vpa_(2.521714853467526e-30) + (-1.366223618221636e-46)+ vpa_(3.141592653589793e-63) + (-1.234567890123457e-79) ...
%     + vpa_(5.678901234567890e-96) + vpa_(9.876543210987654e-112) ) );
b_oct_vpa = vpa_(b_oct);
s_oct = vpa_(a_oct + b_oct);
s_oct_vpa = a_oct_vpa + b_oct_vpa;
err_oct = abs(s_oct - s_oct_vpa);
rel_oct = err_oct / abs(s_oct_vpa);
fprintf('OctoDoubleSlow addition absolute error: %g\n', err_oct);
fprintf('OctoDoubleSlow addition relative error: %g\n', rel_oct);

% Also test OctoDoubleSlow multiplication
p_oct = vpa_(a_oct .* b_oct);
p_oct_vpa = a_oct_vpa .* b_oct_vpa;
err_oct_mul = abs(p_oct - p_oct_vpa);
rel_oct_mul = err_oct_mul / abs(p_oct_vpa);
fprintf('OctoDoubleSlow multiplication absolute error: %g\n', err_oct_mul);
fprintf('OctoDoubleSlow multiplication relative error: %g\n', rel_oct_mul);

