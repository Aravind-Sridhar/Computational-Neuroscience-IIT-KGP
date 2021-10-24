clc;clear all;close all;
tRange=[0 200]; %Range
y0=[1;0];  %Inital conditions

[t,sol45]=ode45(@vdpe1,tRange,y0); %For μ=1 using ODE45
subplot(3,4,1)
plot(t,sol45(:,1));title('ode45 μ=1')
subplot(3,4,2)
plot(sol45(:,1),sol45(:,2));title('ode45 Phase Plane μ=1')

[t1,sol45]=ode45(@vdpe2,tRange,y0); %For μ=0.1 using ODE45
subplot(3,4,3)
plot(t1,sol45(:,1));title('ode45 μ=0.1')
subplot(3,4,4)
plot(sol45(:,1),sol45(:,2));title('ode45 Phase Plane μ=0.1')

[t2,sol45]=ode45(@vdpe3,tRange,y0); %For μ=100 using ODE45
subplot(3,4,5)
plot(t2,sol45(:,1));title('ode45 μ=100')
subplot(3,4,6)
plot(sol45(:,1),sol45(:,2));title('ode45 Phase Plane μ=100')

[t3,sol15s]=ode15s(@vdpe4,tRange,y0); %For μ=1 using ODE15s
subplot(3,4,7)
plot(t3,sol15s(:,1));title('ode15s μ=1')
subplot(3,4,8)
plot(sol15s(:,1),sol15s(:,2));title('ode15s Phase Plane μ=1')

[t4,sol15s]=ode15s(@vdpe5,tRange,y0); %For μ=0.1 using ODE15s
subplot(3,4,9)
plot(t4,sol15s(:,1));title('ode15s μ=0.1')
subplot(3,4,10)
plot(sol15s(:,1),sol15s(:,2));title('ode15s Phase Plane μ=0.1')

[t5,sol15s]=ode15s(@vdpe6,tRange,y0); %For μ=100 using ODE15s
subplot(3,4,11)
plot(t5,sol15s(:,1));title('ode15s μ=100')
subplot(3,4,12)
plot(sol15s(:,1),sol15s(:,2));title('ode15s Phase Plane μ=100')

function dYdt=vdpe1(t,y)
    u=1;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end

function dYdt=vdpe2(t1,y)
    u=0.1;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end

function dYdt=vdpe3(t2,y)
    u=100;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end

function dYdt=vdpe4(t3,y)
    u=100;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end

function dYdt=vdpe5(t4,y)
    u=0.1;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end
function dYdt=vdpe6(t5,y)
    u=1;
    Y1=y(1);
    Y2=y(2);
    dY1dt=Y2;
    dY2dt= u*(1-(Y1^2))*Y2-Y1;
    dYdt=[dY1dt;dY2dt];
end
