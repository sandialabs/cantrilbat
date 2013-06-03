function resid = odeSet(t,x,dxdt)
global A;
resid = dxdt - A*x;