# Coding Standards

(*) Files should have at most one empty line in a row.
(*) All functions and methods & properties blocks should be followed by an empty line. Any block of any kind that contains an empty line should have an empty line before its first indented line, and after its last indented line.
(*) classdef, methods and properties statements should always be followed by an empty line, and have an empty line before their end.
(*) All functions should have an end statement.
(*) All top-level (non-indented) function statements should be followed by an empty line, and they should have an empty line before their end.
(*) Open brackets should be followed by space. Close brackets should be preceeded by space. Colons and other operators should be surrounded by spaces.
(*) Any operations on the left or right side of a colon should be bracketed. I.e. ( a + 1 ) : b NOT a + 1 : b.
(*) Variables should generally have PascalCase. One or two letter lower case variable names are acceptable, but if you want modified versions of them, then use camelCase for the rest of the modified name (e.g. v1Select as a modified version of v1). Underscores in variable names should be avoided unless removing them makes things unclear.
(*) If a function is only called in one place, then it should be deleted and its code should be inserted in that one place.
(*) If a variable is only used in one place, then it should be deleted and its definition should be inserted in that one place.
(*) NaN and Inf should always have that capitalization. (I.e. use NaN, not nan.)
(*) feval or eval should not be used anywhere.
