a=QuadDoubleSlow.rand(200,200);
b=QuadDoubleSlow.rand(200,200);
tic;
c1=a*b;
toc;
a=QuadDouble(a);
b=QuadDouble(b);
tic;
c2=a*b;
toc;
c1=QuadDouble(c1);
disp(max(abs(c1(:)-c2(:))));

