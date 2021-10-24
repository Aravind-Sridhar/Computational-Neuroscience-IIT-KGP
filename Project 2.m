clc;clear;close all;warning off;

Morris_Lecar_MLE

Hodgkin_Huxley_HH

function Morris_Lecar_MLE

global gCa gK gL VCa VK VL phi v1 v2 v3 v4 v5 v6 C Iext;

%Parameters
gCa=4.4; % mS/ cm^2
gK=8; % mS/ cm^2
gL=2; % mS/ cm^2
VCa=120; % mV
VK=-84; % mV
VL=-60; % mV
phi=0.02; 
v1=-1.2; % mV
v2=18; % mV
v3=2; % mV
v4=30; % mV
v5=2; % mV
v6=30;% mV
C=20;  % uF/cm^2
Iext=0; % uA/cm^2


%% Question 2
%Vnull cline and Wnull cline
figure;
hold on
Vnc=@(V)(Iext - gCa *(0.5*(1+tanh((V-v1)/v2)))*(V-VCa)- gL* (V-VL))/(gK* (V-VK)); 
Wnc=@(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc(V)*100, [-80,100],'r');
fplot(@(V) Wnc(V)*100, [-80,100],'b');
xlabel('V(mV)');
ylabel('W');
title('Morris Lecar Phase Plane Plot');

% Equilibrium Points Calculation 

syms V w
Veq= (1/C)*(Iext-gCa*(0.5*(1+tanh((V-v1)/v2)))-gK*w*(V-VK)-gL*(V-VL))==0;
Weq=(0.5*(1+tanh((V-v3)/v4))-w)==0;
Equli_Pt= solve([Veq,Weq],[V,w]);

Veq = double(Equli_Pt.V);
Weq = double(Equli_Pt.w);
grid on
plot(Veq, 100*Weq, 'ko','linewidth', 1.5);

fprintf('Question 2. The equilibrium point of MLE is located at (%d,%d) \n', Veq, Weq);

%Quiver Plot

a=linspace(-80,100,100);
b=linspace(0,0.9,100);

[V_quiv,w_quiv] = meshgrid(a,b);

dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V_quiv-v1)/v2))).*(VCa-V_quiv) + gK*w_quiv .*(VK-V_quiv) + gL*(VL-V_quiv) + Iext);
dw_dt = phi*((0.5*(1+tanh((V_quiv-v3)/v4)))-w_quiv).*cosh((V_quiv-v3)/(2*v4));
quiver(a,b*150,dV_dt,dw_dt*100);
legend('V nullcline', 'w nullcline','Equilibrium Point','Trajectories');

%% Question 3. Stability of Equilibrium Point

syms V w
dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V-v1)/v2)))*(VCa-V) + gK*w*(VK-V) + gL*(VL-V) + Iext);
dw_dt = phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/(2*v4));

J = jacobian([dV_dt, dw_dt],[V,w]);
V = Veq;
w = Weq;
Jmatrix = zeros(2,2);
Jmatrix(1,1) = subs(J(1,1));
Jmatrix(1,2) = subs(J(1,2));
Jmatrix(2,1) = subs(J(2,1));
Jmatrix(2,2) = subs(J(2,2));

eigenValues = eig(Jmatrix);
fprintf('Question 3. The eigen values of the system are %d and %d \n', eigenValues(1), eigenValues(2));
eigenvalue_stability(eigenValues(1),eigenValues(2))

%% Question 5. Action Potential Generation

op = odeset('RelTol',1e-3,'AbsTol',1e-6, 'refine',5, 'MaxStep', 1);

Iext = 0;
tRange = [0, 300];
int = [0,Weq];

phi = 0.01;
[t1, S1] = ode15s(@(t,S)MLE_Solve(t,S),tRange, int, op);

phi = 0.02;
[t2, S2] = ode15s(@(t,S)MLE_Solve(t,S),tRange, int, op);

phi = 0.04;
[t3, S3] = ode15s(@(t,S)MLE_Solve(t,S),tRange, int, op);


figure;
hold on;
plot(t1,S1(:,1));
plot(t2, S2(:,1));
plot(t3, S3(:,1));
xlabel('Time(ms)');
ylabel('Volatage(mV)');
title('Action potentials with different \phi');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

figure;
hold on
Vnc = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc(V), [-80 100],'k');
fplot(@(V) wnc(V), [-80 100],'r');

plot(S1(:,1),S1(:,2));
plot(S2(:,1),S2(:,2));
plot(S3(:,1),S3(:,2));
xlabel('V(mV)');
ylabel('W');
ylim([0,1]);
title('MLE Phase Plane Plot');
legend('V Nullcline','W nullcline','\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

%% Question 6. Depolarization Threshold

Iext = 0;
tSpan = [0, 300];
phi=0.02;
initialV=linspace(-15.2,-14.8,400); 
max_V = zeros(1,400);

flag=1;
for n = 1:400 
    [t,S] = ode15s(@(t,S)MLE_Solve(t,S),tSpan, [initialV(n),Weq], op);
    max_V(n) = max(S(:,1));
    if max_V(n) >= 0 && flag==1
        fprintf("Question 6. Threshold is (%f)",initialV(n));
        threshold = initialV(n);
        flag=0;
    end
end

figure;
hold on
plot(initialV,max_V);
grid on;
xlabel('Initial Voltage(mV)');
ylabel('Maximum Voltage(mV)');
title('Threshold Behavior with change in Initial Voltage');

figure;
hold on
Vnc = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc(V), [-80 100],'m');
fplot(@(V) wnc(V), [-80 100],'k');

V_plot = linspace(threshold-0.1,threshold + 0.1,5);
for i = 1:5
    [t,S] = ode15s(@(t,S)MLE_Solve(t,S),tSpan, [V_plot(i),Weq], op);
    plot(S(:,1),S(:,2));
end
xlabel('V(mV)');
ylabel('W');
ylim([0,1]);
title('MLE Phase Plane Plot for different Initial Voltages around Threshold');
grid on;


%% Question 7. Iext= 86 μA/cm2
Iext = 86;

figure;
hold on
Vnc1 = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc1 = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc1(V), [-80 100]);
fplot(@(V) wnc1(V), [-80 100]);
xlabel('V(in mV)');
ylabel('w');
title('Phase Plane Plot(MLE)');

syms V w
Vnc1_eq = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc1_eq = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
eq_pt1 = solve([Vnc1_eq, wnc1_eq], [V, w]);

V_eq1 = double(eq_pt1.V);
w_eq1 = double(eq_pt1.w);

plot(V_eq1, w_eq1, 'ko', 'linewidth', 2);
text(V_eq1, w_eq1, ['(' num2str(round(V_eq1,5)) ',' num2str(round(w_eq1,5)) ')']);
grid on;

fprintf('\nQuestion 7. The equilibrium point is located at (%d,%d) \n', V_eq1, w_eq1);

dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V-v1)/v2)))*(VCa-V) + gK*w*(VK-V) + gL*(VL-V) + Iext);
dw_dt = phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/(2*v4));

JSym = jacobian([dV_dt, dw_dt],[V,w]);
V = V_eq1;
w = w_eq1;
Jmatrix = zeros(2,2);
Jmatrix(1,1) = subs(JSym(1,1));
Jmatrix(1,2) = subs(JSym(1,2));
Jmatrix(2,1) = subs(JSym(2,1));
Jmatrix(2,2) = subs(JSym(2,2));

eigenValues = eig(Jmatrix);
eigenvalue_stability(eigenValues(1),eigenValues(2))

tSpan = [0,300];
initial = [Veq,Weq];
[t1, S1] = ode15s(@(t,S)MLE_Solve(t,S),tSpan, initial, op);
initial = [V_eq1,w_eq1];
[t2, S2] = ode15s(@(t,S)MLE_Solve(t,S),tSpan, initial, op);
initial = [-27.9,0.17];
[t3, S3] = ode15s(@(t,S)MLE_Solve(t,S),tSpan, initial, op);


plot(S1(:,1),S1(:,2));
plot(S2(:,1),S2(:,2));
plot(S3(:,1),S3(:,2));
ylim([0 1]);
legend('V Nullcline','W Nullcline','Equilibrium point','Eq.Point when Iext = 0','Eq.Point when Iext = 86','Inital Point- Random');

%% Question 8. Running the model backwards in time and True Threshold 
Iext = 86;

figure;
hold on
Vnc1 = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc1 = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc1(V), [-80 100]);
fplot(@(V) wnc1(V), [-80 100]);
grid on
xlabel('V(in mV)');
ylabel('w');
title('Phase Plane Plot(MLE)');


syms V w
Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
eq_pt_1 = solve([Vnc1_eqn, wnc1_eqn], [V, w]);

V_eq1 = double(eq_pt_1.V);
w_eq1 = double(eq_pt_1.w);

plot(V_eq1, w_eq1, 'o');
text(V_eq1, w_eq1, ['(' num2str(round(V_eq1,5)) ',' num2str(round(w_eq1,5)) ')']);
grid on;


tSpan1 = [0,300];
tSpan2 = [0,-1000];

[t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), tSpan1, [Veq, Weq]);
[t2,S2]=ode15s(@(t,S)MLE_Solve(t,S), tSpan1, [-27.9, 0.17]);
[t3,S3]=ode15s(@(t,S)MLE_Solve(t,S), tSpan1, [V_eq1, w_eq1]);
[t4,S4]=ode15s(@(t,S)MLE_Solve(t,S), tSpan2, [-27.9, 0.17]);
plot(S1(:,1), S1(:,2),'g');
plot(S2(:,1), S2(:,2),'y');
plot(S3(:,1), S3(:,2),'b');
plot(S4(:,1), S4(:,2),'m');
legend('W Nullcline','V Nullcline','Equilibrium point', 'Eq.Point for Iext=0','Random Initial point', ...
       'Eq.Point for Iext=86','UPO for Random initial point');

ylim([0 1]);

%% Question 9. Equilibrium Points for Iext = 80, 86, 90 μA/cm2

Iexts = [80, 86, 90];
for i = 1:3
    Iext = Iexts(i);
    fprintf('\nQuestion 9. For Iext = %d \n', Iext);
    % Finding the equilibrium point for this system using MATLAB
    syms V w
    Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    eq_pt_1 = solve([Vnc1_eqn, wnc1_eqn], [V, w]);

    V_eq1 = double(eq_pt_1.V);
    w_eq1 = double(eq_pt_1.w);

    
    fprintf('The Equilibrium point is at (%d,%d)  \n', V_eq1, w_eq1);

    dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V-v1)/v2)))*(VCa-V) + gK*w*(VK-V) + gL*(VL-V) + Iext);
    dw_dt = phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/(2*v4));

    JSymbolic = jacobian([dV_dt, dw_dt],[V,w]);
    V = V_eq1;
    w = w_eq1;
    Jmatrix = zeros(2,2);
    Jmatrix(1,1) = subs(JSymbolic(1,1));
    Jmatrix(1,2) = subs(JSymbolic(1,2));
    Jmatrix(2,1) = subs(JSymbolic(2,1));
    Jmatrix(2,2) = subs(JSymbolic(2,2));

    eigenValues = eig(Jmatrix);
    fprintf('The eigen values are  %f%+fi , %f%+fi \n', real(eigenValues(1)), imag(eigenValues(1)), ...
            real(eigenValues(2)), imag(eigenValues(2)));
    eigenvalue_stability(eigenValues(1),eigenValues(2))
end

figure;
hold on
ylabel('Firing Rate (in Hz)');
xlabel('Iext (uA)');
title('Firing Rate vs External Current');
i = 0;
frequency = zeros(21, 1);
current = zeros(21, 1);
for Iext = 80:1:100
    syms V w
    Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    eq_pt_1 = solve([Vnc1_eqn, wnc1_eqn], [V, w]);
    V_eq1 = double(eq_pt_1.V);
    w_eq1 = double(eq_pt_1.w);
    
    initpoint = [V_eq1 + 0.1, w_eq1 + 0.001];
    tspan = [0 2000];
    [t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), tspan, initpoint);
    i = i+1;
    frequency(i) = freq(t1, S1);
    current(i) = Iext;
    hold on
end
plot(current, frequency);

hold off

%% Question 10. Phase Plane, Null-Clines, Equilibrium Points and Manifolds for new parameters

gCa = 4;
VCa = 120;
gK = 8;
VK = -84;
gL = 2;
VL = -60;
v1 = -1.2;
v2 = 18;
v3 = 12;
v4 = 17.4;
phi = 0.0667;
Iext = 30;
C = 20;


figure;
hold on
Vnc1 = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
wnc1 = @(V) (0.5*(1+tanh((V-v3)/v4)));
fplot(@(V) Vnc1(V), [-80 100], '--');
fplot(@(V) wnc1(V), [-80 100], ':');
xlabel('V(mV)');
ylabel('W');
title('Morris Lecar Phase Plane Plot- New Parameters');


syms V w
Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
soln1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-40, 0]);
soln2 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-20, 0.02]);
soln3 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [0, 0.25]);

V_eq1 = double(soln1.V);
w_eq1 = double(soln1.w);
plot(V_eq1, w_eq1, 'ko');
text(V_eq1, w_eq1, ['(' num2str(round(V_eq1,3)) ',' num2str(round(w_eq1,3)) ')']);
grid on;

V_eq2 = double(soln2.V);
w_eq2 = double(soln2.w);
plot(V_eq2, w_eq2, 'ks');
text(V_eq2, w_eq2, ['(' num2str(round(V_eq2,3)) ',' num2str(round(w_eq2,3)) ')']);
grid on;

V_eq3 = double(soln3.V);
w_eq3 = double(soln3.w);
plot(V_eq3, w_eq3, 'kd');
text(V_eq3, w_eq3, ['(' num2str(round(V_eq3,3)) ',' num2str(round(w_eq3,3)) ')']);
grid on;

fprintf('\nQuestion 10. Equilibrium points are : \n 1. (%f, %f) \n 2. (%f, %f) \n 3. (%f, %f)\n', ...
        V_eq1, w_eq1, V_eq2, w_eq2, V_eq3, w_eq3);

syms V w
V_eq = [V_eq1, V_eq2, V_eq3];
w_eq = [w_eq1, w_eq2, w_eq3];
dV_dt = (1/C)*(gCa*(0.5*(1+tanh((V-v1)/v2)))*(VCa-V) + gK*w*(VK-V) + gL*(VL-V) + Iext);
dw_dt = phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/(2*v4));
JSymbolic = jacobian([dV_dt, dw_dt],[V,w]);
eigVectors = 0;
for i = 1:3
    V = V_eq(i);
    w = w_eq(i);
    Jmatrix = zeros(2,2);
    Jmatrix(1,1) = subs(JSymbolic(1,1));
    Jmatrix(1,2) = subs(JSymbolic(1,2));
    Jmatrix(2,1) = subs(JSymbolic(2,1));
    Jmatrix(2,2) = subs(JSymbolic(2,2));
    eigenValues = eig(Jmatrix);
    if i == 2
        [eigVectors, D] = eig(Jmatrix);
    end
    fprintf('For Equilibrium point %d : The eigen values are  %f%+fi , %f%+fi \n', i, real(eigenValues(1)), imag(eigenValues(1)), ...
            real(eigenValues(2)), imag(eigenValues(2)));
    eigenvalue_stability(eigenValues(1),eigenValues(2))
end


vs = eigVectors(:,1)';
vu = eigVectors(:,2)';


start1 = [V_eq2, w_eq2] +  1 * vs;
start2 = [V_eq2, w_eq2] +  1 * vu;
start3 = [V_eq2, w_eq2] -  1 * vs;
start4 = [V_eq2, w_eq2] -  1 * vu;
[t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), [0 100], start1);
plot(S1(:,1), S1(:,2), 'm');
[t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), [0 -80], start2);
plot(S1(:,1), S1(:,2), 'g');
[t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), [0 100], start3);
plot(S1(:,1), S1(:,2), 'm');
[t1,S1]=ode15s(@(t,S)MLE_Solve(t,S), [0 -30], start4);
plot(S1(:,1), S1(:,2), 'g');
legend('W Nullcline','V Nullcline','Equilibrium Point 1', 'Equilibrium Point 2', 'Equilibrium Point 3', ...
        'Unstable Manifolds', 'Stable Manifolds');

%% Question 11. Rate of firing action potentials


gCa = 4;
VCa = 120;
gK = 8;
VK = -84;
gL = 2;
VL = -60;
v1 = -1.2;
v2 = 18;
v3 = 12;
v4 = 17.4;
phi = 0.0667;
C = 20;

currents = [30, 35, 39, 39.1, 39.2, 39.3, 39.4, 39.5, 39.5, 39.7, 39.8, 39.9, 40, 41, 42,  45,  47, 50];
N = size(currents);
N=N(2);
frequency = zeros(N, 1);
for i = 1:N
    
    Iext = currents(i);
    if (Iext > 37) && (Iext < 42)
        figure;
        hold on
        Vnc1 = @(V) (Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gL*(V-VL))/(gK*(V-VK));
        wnc1 = @(V) (0.5*(1+tanh((V-v3)/v4)));
        fplot(@(V) Vnc1(V), [-80 100], '--');
        fplot(@(V) wnc1(V), [-80 100], ':');
        xlabel('V(in mV)');
        ylabel('w');
        title(strcat('Morris Lecar Phase Plane Plot with Iext = ', num2str(Iext)));
        hold on;
    end
    
    syms V w
    Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    soln1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-40, 0]);
    soln2 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-20, 0]);
    soln3 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [0, 0.2]);

    V_eq1 = double(soln1.V);
    w_eq1 = double(soln1.w);
    V_eq2 = double(soln2.V);
    w_eq2 = double(soln2.w);
    V_eq3 = double(soln3.V);
    w_eq3 = double(soln3.w);
    
    if isequal([V_eq1, w_eq1], [V_eq2, w_eq2]) 
        fprintf("\nQuestion 11. For Iext = %f, there is only one distinct equilibrium point\n", Iext);

        [t,S]=ode15s(@(t,S)MLE_Solve(t,S), [0 3000], [V_eq3 + 0.1, w_eq3+0.01]);
        if (Iext > 37) && (Iext < 42)
            plot(V_eq3, w_eq3, 'o');
            text(V_eq3, w_eq3, ['(' num2str(round(V_eq3,3)) ',' num2str(round(w_eq3,3)) ')']);
            grid on;
            plot(S(:,1), S(:,2));
        end
        frequency(i) = freq(t, S);
    else
        fprintf("\nQuestion 11. For Iext = %f, there are three distinct equilibrium points\n", Iext);
        [t,S]=ode15s(@(t,S)MLE_Solve(t,S), [0 3000], [V_eq2 - 0.1, w_eq2-0.01]);
        frequency(i) = 0;
        
        if (Iext >= 37) && (Iext < 42)
            plot(V_eq1, w_eq1, 'o');
            grid on;
            plot(V_eq2, w_eq2, '-k');
            grid on;
            plot(V_eq3, w_eq3, 'o');
            text(V_eq3, w_eq3, ['(' num2str(round(V_eq3,3)) ',' num2str(round(w_eq3,3)) ')']);
            grid on;
            plot(S(:,1), S(:,2));
        end
    end

end
hold off

figure;
hold on
ylabel('Firing Rate (Hz)');
xlabel('Iext (uA)');
title('Firing Rate vs External Current');
plot(currents, frequency);
hold off    

end


function dS = MLE_Solve(t,S)

global C;
global gCa;
global VCa;
global gK;
global VK;
global gL;
global VL;
global v1;
global v2;
global v3;
global v4;
global phi;
global Iext;

%Defining State Variables
V=S(1);
w=S(2);

%Local functions
m_inf = (0.5*(1+tanh((V-v1)/v2)));
w_inf = (0.5*(1+tanh((V-v3)/v4)));

ddt_V = (1/C)*(gCa*m_inf*(VCa-V) + gK*w*(VK-V) + gL*(VL-V)+Iext);
ddt_w = phi*(w_inf-w)*cosh((V-v3)/(2*v4));

dS=[ddt_V; ddt_w];

end
function F = freq(t, S)
    n = size(S);
    n = n(1);
    S = S(:, 1);
    % Spike generation check:
    has_spikes = 0;
    for i = 1:n
        if S(i) > 0
            has_spikes = 1;
        end
    end
    if has_spikes == 0
        F = 0;
        return
    end
    
    
    while S(n) > 0 
        n = n - 1;
    end
    
    while S(n) < 0 
        n = n -1 ;
    end
    t2 = t(n);
    
    while S(n) >0
        n = n-1;
    end
   
    while S(n) < 0
        n = n-1;
    end
    t1 = t(n);
    F = 1000 / (t2 - t1);
end
function eigenvalue_stability(v1,v2)

if isreal(v1)
    if v1*v2>0
        if v1<0
        display('The point is stable')
        end
        if v1>0
        display('The point is stable')
        end
    end
    if v1*v2<0
        display('The point is a saddle node')
    end
else 
    if real(v1)>0
        display('The point is unstable and diverging spiral')
    end
if real(v1)<0
    display('The point is stable and converging spiral')
end
end
end

function Hodgkin_Huxley_HH

%Question 12. Hodgkin Huxley Model

    global C Iext gK gNa gL VK VNa VL eps f hconst;
    
    C = 1;
    gK = 36;
    VK = -72;
    gNa = 120;
    VNa = 55;
    gL = 0.3;
    Iext = 0;
    eps = 1e-10;
    
    options = odeset('RelTol',1e-9,'AbsTol',1e-9, 'refine',5, 'MaxStep', 1);
    %% Question 13. For Vr = -60 
   
    
    Vr = -60;
    alphan = -0.01 * (Vr + eps + 50)/(exp(-(Vr + eps + 50)/10)-1);
    alpham = -0.1 * (Vr + eps + 35)/(exp(-(Vr + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(Vr + 60)/20);
    betan = 0.125 * exp(-(Vr + 60)/80);
    betam = 4 * exp(-(Vr + 60)/18);
    betah = 1/(exp(-(Vr + 30)/10) + 1);
    
    mInf = alpham/(alpham + betam);
    nInf = alphan/(alphan + betan);
    hInf = alphah/(alphah + betah);
    hconst = hInf;
 
    E_leak = Vr - (1/gL)*(Iext - gK * (nInf^4) * (Vr - VK) - gNa * (mInf^3) * hInf * (Vr - VNa));
    
    fprintf('\nQuestion 13. EL = %d \n', E_leak);
    
    VL = E_leak;
    
    %For Iext=10uA
    Iext = 10;
    tSpan = [0 100];
    
    Sinit = [-60, nInf, mInf, hInf]; 
    [t1, S] = ode15s(@(t,S)HH(t,S), tSpan, Sinit, options);
    figure;
    plot(t1, S(:,1));
    xlabel('Time(ms)');
    ylabel('Voltage (mV)');
    title('HH model with given parameters');
    
    %% Question 14. Stability and Threshold
    
    Iext = 0;
    impulseI = linspace(0,15,100);
    for i = 1:100 
        Sinit = [-60+impulseI(i)/C, nInf, mInf, hInf];
        [t1, S] = ode15s(@(t,S)HH(t,S), tSpan, Sinit, options);
        max_V(i) = max(S(:,1));
    end
    thres = (min(max_V(:)) + max(max_V(:)))/2;
    thres_i = 500;
    for i = 1:100
        if max_V(i) >= thres && i < thres_i
            thres_i = i;
        end
    end
    fprintf("\nQuestion 14. Threshold: %f mA\n",impulseI(thres_i));
    figure;
    plot(impulseI, max_V);
    grid on;
    xlabel('Impulse current(mA)');
    ylabel('Peak Voltage (mV)');
    title('Current pulse Threshold');
    
    fprintf('Question 14. Iext =0\n');
    stability(0);
    
%% Question 15. For Iext between 8uA and 12uA

    for i=8:12
        fprintf("\nQuestion 15. Iext = %d \n", i);
        stability(i);
    end
 %% Question 16. Myotonia (Na+ channel inactivation)
figure;
Sinit = [-52, nInf, mInf, hInf];
for f = [0, 0.1, 0.17, 0.2]
     [t1, S] = ode15s(@(t,S)HH_Inac(t,S), tSpan, Sinit, options);
     hold on;
      plot(t1, S(:,1));
 end
xlabel('Time(ms)');
ylabel('Voltage (mV)');
title('Action Potential in the case of Myotonia due to Na+ channel inactivation');
legend('f = 0','f = 0.1','f = 0.17','f = 0.2');
grid on;

%% Question 17. V-n reduced form of HH model

f = 0;
Iext = 10;
tSpan = [0 100];
Sinit = [-60, nInf]; 
[t1, S] = ode15s(@(t,S)Reduced_HH(t,S), tSpan, Sinit, options);
figure;
plot(t1, S(:,1));
xlabel('Time(ms)');
ylabel('Voltage (mV)');
title('V-n Reduced HH model-Current injection');
hold off

Iext = 0;
f = 0;
impulseI = linspace(0,15,5);
figure;
xlabel('Time(ms)');
ylabel('Voltage (mV)');
title('V-n Reduced HH model-Current pulse');
hold on
for i = 1:5
Sinit = [-60+impulseI(i)/C, nInf];
[t1, S] = ode15s(@(t,S)Reduced_HH(t,S), [0 20], Sinit, options);
plot(t1, S(:,1));
hold on
end
legend(strcat('Iimpulse = ', num2str(impulseI(1))), strcat('Iimpulse = ', num2str(impulseI(2))), ...
strcat('Iimpulse = ', num2str(impulseI(3))), strcat('Iimpulse = ', num2str(impulseI(4))), ...
strcat('Iimpulse = ', num2str(impulseI(5))));
hold off

%% Question 18. Phase Plane analysis of V-n reduced HH model
fs = linspace(0.02, 0.4, 7);
    Iext = 0;
 
for f = fs
figure;
hold on
syms V
alphn = @(V) -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
alphm = @(V) -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
betn = @(V) 0.125 * exp(-(V + 60)/80);
betm = @(V) 4 * exp(-(V + 60)/18);
minf = @(V) alphm(V) / (alphm(V) + betm(V));
    
Vnc = @(V) ( (Iext - gNa*(1-f)*(minf(V)^3)*hconst*(V-VNa) - gNa*f*(minf(V)^3)*(V-VNa) - gL*(V-VL))/(gK * (V-VK)) )^(1/4);
nnc = @(V) alphn(V)/(alphn(V) + betn(V));                                                                   
fplot(@(V) Vnc(V), [-80 100], ':','linewidth',1.2);
fplot(@(V) nnc(V), [-80 100], '-.');
xlabel('V(mV)');
ylabel('n');
title(strcat('Phase Plane Plot(V-n Reduced HH model) with f = ', num2str(f)));
hold on;
        
impulseI = linspace(0, 15, 5);
for i=1:5
Sinit = [-60+impulseI(i)/C, nInf];
[t1, S1] = ode15s(@(t,S)Reduced_HH(t,S), [0, 500], Sinit, options);
hold on;
plot(S1(:,1), S1(:,2));
end
legend('V-null cline', 'n-null cline', 'Impulse=0' , 'Impulse=3.75', ...
                'Impulse=7.5', 'Impulse=11.25', 'Impulse=15');
end

%% Question 19. Anode Break Excitation
figure;
tSpan1 = [0 20];
Sinit1 = [-60, nInf, mInf, hInf];
Iext = -3;
[t1, S1] = ode15s(@(t,S)HH(t,S), tSpan1, Sinit1, options);
tSpan2 = [20 100];
Sinit2 = [S1(end,1), S1(end,2), S1(end,3), S1(end,4)];
nAfterA = S1(end,2);
hAfterA = S1(end,4);
Iext = 0;
[t2, S2] = ode15s(@(t,S)HH(t,S), tSpan2, Sinit2, options);
totalT = [t1; t2];
totalV = [S1(:,1); S2(:,1)];
plot(totalT, totalV);
xlabel('Time(ms)');
ylabel('Voltage(mV)');
title('Action potential in case of Anode Break excitation');
grid on;
    
%% Question 20. Anode break of V-m reduced HH model 
    
Vr = -60;
    
alphan = @(V) -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
alpham = @(V) -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
alphah = @(V) 0.07 * exp(-(V + 60)/20);
betan = @(V) 0.125 * exp(-(V + 60)/80);
betam = @(V) 4 * exp(-(V + 60)/18);
betah = @(V) 1/(exp(-(V + 30)/10) + 1);
    
mInf = alpham(Vr)/(alpham(Vr) + betam(Vr));
nInf = alphan(Vr)/(alphan(Vr) + betan(Vr));
hInf = alphah(Vr)/(alphah(Vr) + betah(Vr));
    
Iext = -3;
tSpan1 = [0 20];
Sinit1 = [-60, mInf];
[t1, S1] = ode15s(@(t,S)HH_Vm_Reduced(t,S,nInf,hInf), tSpan1, Sinit1, options);
Vnc1 = @(V) (((Iext - gK * (V -VK) * (nInf^4)- gL*(V-VL))/(gNa*hInf*(V-VNa)))^(1/3));
mnc = @(V) alpham(V)/(alpham(V) + betam(V));
figure;
hold on;
fplot(@(V) Vnc1(V), [-80 100], ':');
fplot(@(V) mnc(V), [-80 100], '--');
    
Iext = 0;
tSpan2 = [20 100];
Sinit2 = [S1(end,1), S1(end,2)];
    
nInf1 = nAfterA;
hInf1 = hAfterA;
    
[t2, S2] = ode15s(@(t,S)HH_Vm_Reduced(t,S,nInf1,hInf1), tSpan2, Sinit2, options);
Vnc2 = @(V) ( (Iext - gK * (V -VK) * (nInf^4)- gL*(V-VL))/(gNa*hInf*(V-VNa)))^(1/3);
fplot(@(V) Vnc2(V), [-80 100], '-*');
    
totalV = [S1(:,1); S2(:,1)];
totalm = [S1(:,2); S2(:,2)];                                                              
plot(totalV,totalm);
ylabel('m');
xlabel('Voltage(mV)');
title('Anode Break Phase Plane Plot');
grid on;
  
Iext = -3;
syms V m;
    
V_nc = (1/C) * (Iext - gK * nInf^4 * (V - VK) - gNa * m^3 * hInf* (V - VNa) - gL * (V - VL)) == 0;
m_nc = alpham *(1-m) - betam*m == 0 ;
    
eq_pt = solve([V_nc, m_nc], [V, m]);
V_eq1 = double(eq_pt.V);
m_eq1 = double(eq_pt.m);
fprintf('\nQuestion 20. Case 1-Equilibrium Point V=%f m=%f\n',eq_pt.V,eq_pt.m);
    
syms V m;
    
Iext=0;
V_nc = (1/C) * (Iext - gK * nInf1^4 * (V - VK) - gNa * m^3 * hInf1* (V - VNa) - gL * (V - VL)) == 0;
m_nc = alpham *(1-m) - betam*m == 0 ;
    
eq_pt = solve([V_nc, m_nc], [V, m]);
V_eq1 = double(eq_pt.V);
m_eq1 = double(eq_pt.m);
    
fprintf('Question 20. Case 2-Equilibrium Point V=%f m=%f\n',eq_pt.V,eq_pt.m);

end

function dS = HH(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * m^3 * h * (V - VNa) - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    ddt_m = alpham * (1 - m) - betam * m;
    ddt_h = alphah * (1 - h) - betah * h;
    
    dS = [ddt_V; ddt_n; ddt_m; ddt_h];
end

function stability(Iext)
    global C;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    
    syms V n m h 
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    mInf = alpham/(alpham + betam);
    nInf = alphan/(alphan + betan);
    hInf = alphah/(alphah + betah);
    
    V_nc = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * m^3 * h* (V - VNa) - gL * (V - VL)) == 0;
    n_nc = alphan *(1-n) - betan*n == 0 ;
    m_nc = alpham * (1 - m) - betam * m ==0 ;
    h_nc = alphah * (1 - h) - betah * h ==0 ;
    
    eq_pt = solve([V_nc, n_nc, m_nc, h_nc], [V, n, m ,h]);
    V_eq1 = double(eq_pt.V);
    n_eq1 = double(eq_pt.n);
    m_eq1 = double(eq_pt.m);
    h_eq1 = double(eq_pt.h);
    
    for p=1:length(eq_pt)
        fprintf('Equilibrium Point V=%f n=%f m=%f h=%f\n',eq_pt(p).V,eq_pt(p).n,eq_pt(p).m,eq_pt(p).h);
    end
 
    dV_dt = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * m^3 * h* (V - VNa) - gL * (V - VL));
    dn_dt = alphan *(1-n) - betan*n;
    dm_dt = alpham * (1 - m) - betam * m ;
    dh_dt = alphah * (1 - h) - betah * h ;
    
    JSymbolic = jacobian([dV_dt, dn_dt, dm_dt, dh_dt],[V,n,m,h]);
    V = V_eq1;
    n = n_eq1;
    m = m_eq1;
    h = h_eq1;
    Jmatrix = zeros(4,4);
    Jmatrix(1,1) = subs(JSymbolic(1,1));
    Jmatrix(1,2) = subs(JSymbolic(1,2));
    Jmatrix(1,3) = subs(JSymbolic(1,3));
    Jmatrix(1,4) = subs(JSymbolic(1,4));
    Jmatrix(2,1) = subs(JSymbolic(2,1));
    Jmatrix(2,2) = subs(JSymbolic(2,2));
    Jmatrix(2,3) = subs(JSymbolic(2,3));
    Jmatrix(2,4) = subs(JSymbolic(2,4));
    Jmatrix(3,1) = subs(JSymbolic(3,1));
    Jmatrix(3,2) = subs(JSymbolic(3,2));
    Jmatrix(3,3) = subs(JSymbolic(3,3));
    Jmatrix(3,4) = subs(JSymbolic(3,4));
    Jmatrix(4,1) = subs(JSymbolic(4,1));
    Jmatrix(4,2) = subs(JSymbolic(4,2));
    Jmatrix(4,3) = subs(JSymbolic(4,3));
    Jmatrix(4,4) = subs(JSymbolic(4,4));
    
    eigenValues = eig(Jmatrix);
    fprintf('The eigen values are %f%+fi , %f%+fi , %f%+fi , %f%+fi \n', real(eigenValues(1)), imag(eigenValues(1)), ...
            real(eigenValues(2)), imag(eigenValues(2)), real(eigenValues(3)), imag(eigenValues(3)), ...
            real(eigenValues(4)), imag(eigenValues(4)));
    k=0;
    for x=1:4
        if real(eigenValues(x)) < 0 
            k = k+1;
        else
            k = k-1;
        end
    end
    if k == 4 
        fprintf("Stable\n");
    elseif k == -4
        fprintf("Unstable\n");
    else
        fprintf("Stability cannot be determined as 4-D plot is required\n");
    end
    
end
function dS = HH_Inac(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    global f;
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * (1-f) * m^3 * h * (V-VNa) - gNa * f * m^3 * (V-VNa)  - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    ddt_m = alpham * (1 - m) - betam * m;
    ddt_h = alphah * (1 - h) - betah * h;
    
    dS = [ddt_V; ddt_n; ddt_m; ddt_h];
end
function dS = Reduced_HH(t,S)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    global f;
    global hconst;
    V = S(1);
    n = S(2);
    
    alphan =  -0.01 * (V + eps + 50)/(exp(-(V + eps + 50)/10)-1);
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    alphah = 0.07 * exp(-(V + 60)/20);
    betan = 0.125 * exp(-(V + 60)/80);
    betam = 4 * exp(-(V + 60)/18);
    betah = 1/(exp(-(V + 30)/10) + 1);
    
    mInf = alpham/(alpham + betam);
    
    ddt_V = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * (1-f) * mInf^3 * hconst * (V-VNa) - gNa * f * mInf^3 * (V-VNa)  - gL * (V - VL));
    ddt_n = alphan * (1 - n) - betan * n;
    
    dS = [ddt_V; ddt_n];
end
function dS = HH_Vm_Reduced(t,S,nInf,hInf)
    global C;
    global Iext;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global eps;
    
    V = S(1);
    m = S(2);
    
    alpham =  -0.1 * (V + eps + 35)/(exp(-(V + eps + 35)/10)-1);
    betam = 4 * exp(-(V + 60)/18);
    
    ddt_V = (1/C) * (Iext - gK * nInf^4 * (V - VK) - gNa * m^3 * hInf * (V - VNa) - gL * (V - VL));
    ddt_m = alpham * (1 - m) - betam * m;
    
    dS = [ddt_V; ddt_m];
end