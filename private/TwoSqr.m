function [ p1, p2 ] = TwoSqr( a )
    p1 = a .* a;
    [ a1, a2 ] = Split( a );
    p2 = ( ( a1 .* a1 - p1 ) + 2 .* a1 .* a2 ) + a2 .* a2;
end
