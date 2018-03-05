function p = fitShiftedExp(x,y)



ft1 = fittype( 'shiftedExpCurveModel(x,a,b,c,d)' );
ftopt1 = fitoptions(ft1);
ftopt1.Lower=[0, 0, min(x), 0];
ftopt1.Upper=[Inf Inf 0 0.3];
% max(y((end-10):end))
f1 = fit( x, y, ft1,ftopt1);
figure
plot( f1, x, y ) ;

p=[f1.a f1.b f1.c f1.d];


end