clear all
global A;
eps = 1e-1;
A = [-100 99/eps; 0 -eps]./2;
title{1} = 'Experiments of Sharma et al. 1977';
title{2} = 'Cantera Simulation, R_{Diffusive} = 0.04 Ohm';
files = dir('outfile*');
figure(1)
% hold on
i = 1;
    data = textread(files(i).name,'','delimiter',',','headerlines',1);
    time = data(:,1);
    value1 = data(:,2);
    value2 = data(:,3);
    %     voltage = 0.194+voltage-voltage(1);
    %voltage = 0.04*current(i)+voltage;
    h(i) = plot(time,value1,'b*');
    %    hh(i) = plot(capacity/1.825,voltage,'--','Color',colors(i,:));
    hold on
    hh(i) = plot(time,value2,'r*');

tspan = [0:0.01:1];

[v d] = eig(A);
A*v(:,1);
y1 = v(1,1)*exp(d(1,1)*tspan) + v(1,2)*exp(d(2,2)*tspan);
y2 = v(2,1)*exp(d(1,1)*tspan) + v(2,2)*exp(d(2,2)*tspan);
% plot(tspan, y1, 'b')
% plot(tspan, y2, 'r')


yinit = [y1(1); y2(1)]
yinit = [0; 10*eps]
yprimeinit = [0; 0];
[tout yout] = ode15i(@odeSet,tspan,yinit,yprimeinit);

plot(tspan, yout(:,1), 'bo')
plot(tspan, yout(:,2), 'ro')


hold off