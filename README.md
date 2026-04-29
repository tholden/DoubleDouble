# Extended Precision Numerical Library

A comprehensive MATLAB library for extended precision arithmetic, providing progressively higher precision tiers starting from double-double (~32 digits) up to oct-double (~128 digits) precision.

This library was originally inspired by the QD library (<http://crd-legacy.lbl.gov/~dhbailey/mpdist/>), but has been modernized and expanded into a full polymorphic class hierarchy (`BaseExtDouble`) with native MATLAB integration, optimized recursive arithmetic, full array support, and experimental complex number capabilities.

## The Numerical Hierarchy

The library provides four extended precision classes, each building upon the last to achieve higher precision bounds:

1. **`DoubleDouble`**: The foundational class. Represents numbers as an unevaluated sum of two 53-bit `double` primitives, providing roughly 106 bits (~32 decimal digits) of precision.
2. **`QuadDouble`**: The standard quad-double implementation. Provides 212 bits (~64 decimal digits) of precision using explicitly overloaded, optimized numerical kernels.
3. **`QuadDoubleSlow`**: A mathematically equivalent pure-MATLAB fallback to `QuadDouble`. It achieves 212-bit precision by recursively executing Dekker's split/TwoProd algorithms using `DoubleDouble` objects as its underlying components rather than explicit scalar kernels.
4. **`OctDouble`**: The highest precision class available in the hierarchy. Provides 424 bits (~128 decimal digits) of precision by recursively leveraging `QuadDouble` objects as its fundamental underlying representation.

## Features

The classes natively support almost all standard MATLAB operations, seamlessly integrating into existing codebases. Because they inherit from `BaseExtDouble`, operations scale polymorphicly across all precision tiers:

- **Basic Arithmetic**: `+`, `-`, `.*`, `./`, `.\`, `.^`
- **Linear Algebra**: `lu`, `qr`, `chol`, `ldl`, `inv`, `det`, `eig`, `mldivide` (`\`), `mrdivide` (`/`), `dot`, `norm`
- **Trigonometry & Exponentials**: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `exp`, `log`, `log2`, `log10`
- **Statistics & Arrays**: `sum`, `prod`, `mean`, `median`, `std`, `var`, `cumsum`, `diff`, `reshape`, `repmat`, `sort`, `unique`
- **Relational Operations**: `<`, `>`, `<=`, `>=`, `==`, `~=`, `any`, `all`, `find`
- **Complex Support**: `real`, `imag`, `conj`, `angle`

## Testing & Validation

The library is backed by a rigorous unit testing suite located in the `UnitTests/` directory.

To execute the test suite:

```matlab
cd UnitTests/
res = RunTests()
```

The suite thoroughly validates standard arithmetic, matrix decompositions, and transcendental approximations. A critical component is the `TestCrossValidationVPA` test suite, which isolates rounding errors by evaluating exact binary conversions (`vpa(sym(Data, 'f'), 135)`) against MATLAB's Symbolic Math Toolbox (135-digit ground truth). This enforces strict error bounds to guarantee precision boundaries across the hierarchy:

- `DoubleDouble`: Enforces strict accuracy to a `1e-30` absolute tolerance.
- `QuadDouble` & `QuadDoubleSlow`: Enforces strict accuracy to a `1e-60` absolute tolerance.
- `OctDouble`: Enforces strict accuracy to a `1e-120` absolute tolerance.

## Development & Coding Standards

If you are extending or modifying the library, you should strictly adhere to the project conventions outlined in [CODING_STANDARDS.md](CODING_STANDARDS.md).

These standards enforce:

- Architectural encapsulation (e.g., using `.Assign()` and `.Index()` rather than direct native subsetting).
- Strict casing conventions (PascalCase for files, static initializers, and methods).
- Mandatory spacing, bracket padding, and blank-line declarations.

All new logic must successfully clear the `RunTests` execution with zero regressions before being considered mathematically sound.
