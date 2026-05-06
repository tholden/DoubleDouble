function [ S1, S2, S3 ] = DDNormalize( A1, A2, A3 )
    S1 = A1 + A2;
    T = S1 - A1;
    S2 = A2 - T;
    if nargin >= 3
        T2 = S2 + A3;
        T3 = T2 - S2;
        S3 = A3 - T3;
        S2 = T2;
    end
end
