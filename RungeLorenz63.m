function [K4] = RungeLorenz63(X,tstep,a)

K1 = tstep*Lorenz63Dx(X,a);
K2 = tstep*Lorenz63Dx(X+.5*K1,a);
K3 = tstep*Lorenz63Dx(X+.5*K2,a);
K4 = tstep*Lorenz63Dx(X+K3,a);

K4 = (1/6)*(K1+2*K2+2*K3+K4);

return;


