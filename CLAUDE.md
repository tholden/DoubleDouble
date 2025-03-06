# DoubleDouble Project Guidelines

## Testing Commands
- Run all tests: `Results = RunTests()`
- Run a single test: `run( DoubleDoubleTest( 'TestMethodName' ) )`
  Example: `run( DoubleDoubleTest( 'TestConstructorScalar' ) )`

## Coding Standards
- **Spacing**: Space after opening brackets and before closing brackets
  - Correct: `[ 1, 2, 3 ]` and `function Test( Arg )`
  - Incorrect: `[1, 2, 3]` or `function Test(Arg)`
- **Operators**: Always surrounded by spaces: `A + B`, not `A+B`
- **Array indexing**: Always use spaces: `A( 1 )`, not `A(1)`
  - Exception: Use `(:)` without spaces
- **Naming**:
  - Variables: PascalCase for vars longer than 2 letters
  - Functions: PascalCase (TestConstructorEmpty, not testConstructorEmpty)

## Function Structure
- Include header comments with function purpose and parameters
- Follow MATLAB's OOP syntax: `ToSumOfDoubles( A )` not `A.ToSumOfDoubles()`
- Use appropriate tolerances for numerical comparisons: `AbsTol = 1e-30` and `RelTol = 1e-15`
