DoubleDouble
============
A MATLAB library for extended ("double double") precision, giving close to quad precision.

Based on the QD library available from: http://crd-legacy.lbl.gov/~dhbailey/mpdist/

The DoubleDouble class has the following public properties:
 * eps
 * pi
 
The DoubleDouble class has the following public, non-static methods:
 * DoubleDouble
 * disp
 * ToSumOfDoubles
 * double
 * isreal
 * isfinite
 * isinf
 * isnan
 * real
 * imag
 * conj
 * angle
 * size
 * length
 * numel
 * numArgumentsFromSubscript
 * end
 * repmat
 * reshape
 * isequal
 * isempty
 * diag
 * tril
 * triu
 * plus
 * minus
 * uminus
 * uplus
 * times
 * mtimes
 * rdivide
 * ldivide
 * mldivide
 * mrdivide
 * power
 * mpower
 * lt
 * gt
 * le
 * ge
 * ne
 * eq
 * colon
 * ctranspose
 * transpose
 * horzcat
 * vertcat
 * subsref
 * subsasgn
 * subsindex
 * sum
 * prod
 * max
 * min
 * cumsum
 * diff
 * cumprod
 * cummax
 * cummin
 * dot
 * norm
 * abs
 * sign
 * floor
 * ceil
 * fix
 * round
 * realsqrt
 * sqrt
 * sqrtm
 * exp
 * expm
 * log
 * log2
 * log10
 * logm
 * funm
 * sincos
 * sin
 * asin
 * cos
 * acos
 * tan
 * atan
 * atan2
 * sinh
 * asinh
 * cosh
 * acosh
 * sinhcosh
 * tanh
 * atanh
 * mod
 * rem
 * lu
 * qr
 * det
 * inv
 * chol
 * ldl
 * eig

The DoubleDouble class has the following public, static methods:
 * IsEqualWithExpansion
 * ones
 * zeros
 * eye
 * nan
 * inf
 * rand
 * randn
 * randi
 * Plus
 * Minus
 * Times
 * MTimes
 * RDivide
 * LDivide
 * MLDivide
 * MRDivide
 * Sum
 * CumSum
 * Diff
 * Prod
 * CumProd
 * Norm
 * Max
 * CumMax
 * Min
 * CumMin
 * Dot
 * ExpandSingleton
