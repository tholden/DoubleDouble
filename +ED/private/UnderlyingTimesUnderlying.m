function [ p1, p2 ] = UnderlyingTimesUnderlying( a, b )
    if isa( a, 'ED.BaseQuadDouble' )
        a1_obj = DoubleDouble(a.v1); a2_obj = DoubleDouble(a.v2);
        a0 = a1_obj.v1; a1 = a1_obj.v2; a2 = a2_obj.v1; a3 = a2_obj.v2;
        b1_obj = DoubleDouble(b.v1); b2_obj = DoubleDouble(b.v2);
        b0 = b1_obj.v1; b1 = b1_obj.v2; b2 = b2_obj.v1; b3 = b2_obj.v2;

        qa0 = a.Promote(a0); qa1 = a.Promote(a1); qa2 = a.Promote(a2); qa3 = a.Promote(a3);
        qb0 = b.Promote(b0); qb1 = b.Promote(b1); qb2 = b.Promote(b2); qb3 = b.Promote(b3);

        p00 = qa0 .* qb0; p01 = qa0 .* qb1; p10 = qa1 .* qb0; p11 = qa1 .* qb1;
        p02 = qa0 .* qb2; p20 = qa2 .* qb0; p12 = qa1 .* qb2; p21 = qa2 .* qb1;
        p03 = qa0 .* qb3; p30 = qa3 .* qb0; p13 = qa1 .* qb3; p31 = qa3 .* qb1;
        p22 = qa2 .* qb2; p23 = qa2 .* qb3; p32 = qa3 .* qb2; p33 = qa3 .* qb3;
        
        z = a.Promote(0);

        [s1, e1] = DDPlusDD(p00, z, p01, z);
        [s2, e2] = DDPlusDD(p10, z, p11, z);
        [s3, e3] = DDPlusDD(s1, e1, s2, e2);

        [s4, e4] = DDPlusDD(p02, z, p20, z);
        [s5, e5] = DDPlusDD(p12, z, p21, z);
        [s6, e6] = DDPlusDD(s4, e4, s5, e5);

        [s7, e7] = DDPlusDD(s3, e3, s6, e6);

        [s8, e8] = DDPlusDD(p03, z, p30, z);
        [s9, e9] = DDPlusDD(p13, z, p31, z);
        [s10, e10] = DDPlusDD(s8, e8, s9, e9);

        [s11, e11] = DDPlusDD(p22, z, p23, z);
        [s12, e12] = DDPlusDD(p32, z, p33, z);
        [s13, e13] = DDPlusDD(s11, e11, s12, e12);

        [s14, e14] = DDPlusDD(s10, e10, s13, e13);
        [p1, p2] = DDPlusDD(s7, e7, s14, e14);
        return
    end

    p1 = a .* b;
    [ a1, a2 ] = Split( a );
    [ b1, b2 ] = Split( b );
    [ t1, e1 ] = UnderlyingPlusUnderlying( a1 .* b1, -p1 );
    [ t2, e2 ] = UnderlyingPlusUnderlying( t1, a1 .* b2 );
    [ t3, e3 ] = UnderlyingPlusUnderlying( t2, a2 .* b1 );
    p2 = t3 + (a2 .* b2 + e1 + e2 + e3);
end
