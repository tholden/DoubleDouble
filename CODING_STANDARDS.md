# Coding Standards

* There should be at most one statement per line.
* All functions should have an end statement.
* Files should have at most one empty line in a row.
* Multi-line comments should be surrounded by empty lines.
* All functions and methods & properties blocks should be followed by an empty line.
* Any block of any kind that contains an empty line should have an empty line before its first indented line, and after its last indented line.
* Indenting should use steps of 4 spaces.
* Any block in which the contents are indented by either one or two steps should have an empty line before the first indented line and after the last indented line. (For example, any classdef, methods and properties statements should always be followed by an empty line, and have an empty line before their end.)
* Open brackets should be followed by space. Close brackets should be preceeded by space. Colons and other operators should be surrounded by spaces.
* Any operations on the left or right side of a colon should be bracketed. I.e. ( a + 1 ) : b NOT a + 1 : b.
* Variables should generally have PascalCase. One or two letter lower case variable names are acceptable, but if you want modified versions of them, then use camelCase for the rest of the modified name (e.g. v1Select as a modified version of v1). Underscores in variable names should be avoided unless removing them makes things unclear.
* NaN and Inf should always have that capitalization. (I.e. use NaN, not nan.)
* feval or eval should not be used anywhere.
