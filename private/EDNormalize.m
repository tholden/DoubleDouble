function [ s1, s2 ] = EDNormalize( a1, a2 )
    s1 = a1 + a2;
    t = s1 - a1;
    s2 = a2 - t;
end
