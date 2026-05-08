function [ s1, s2 ] = UnderlyingPlusUnderlying( a, b ) % AKA two_sum
    a = Normalize( a );
    b = Normalize( b );
    s1 = a + b;
    bb = s1 - a;
    t11 = s1 - bb;
    t2 = b - bb;
    t1 = a - t11;
    s2 = t1 + t2;
end
