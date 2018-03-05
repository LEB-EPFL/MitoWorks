function p = fitSigmoid(x,y)



ft1 = fittype( 'sigmoidCurveModel(x,a,b,c,d)' );
ftopt1 = fitoptions(ft1);
ftopt1.Lower=[0, 0.7, min(x), 0];
ftopt1.Upper=[0.3 1 0 Inf];
% max(y((end-10):end))
f1 = fit( x, y, ft1,ftopt1);
figure
plot( f1, x, y ) ;

p=[f1.a f1.b f1.c f1.d];

ft2 = fittype( 'fixedSigmoidCurveModel(x,a,b,c)' );
ftopt2 = fitoptions(ft2);
ftopt2.Lower=[f1.a,f1.b, min(x)];
ftopt2.Upper=[f1.a f1.b 0];
f2 = fit( x, y, ft2,ftopt2);
figure
plot( f2, x, y ) ;

d=log10(((f2.b-f2.a)/(0.99*f2.b-f2.a))-1)/f2.c;

ft3 = fittype( 'sigmoidCurveModel(x,a,b,c,d)' );
ftopt3 = fitoptions(ft3);
ftopt3.Lower=[f1.a,f1.b, f2.c d];
ftopt3.Upper=[f1.a f1.b f2.c d];
f3 = fit( x, y, ft3,ftopt3);
figure
% plot( f3, x, y ) ;

end