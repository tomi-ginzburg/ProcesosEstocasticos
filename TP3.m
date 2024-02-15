clear all; close all; clc;

%% EJERCICIO 1

N = 20000;
m = 500; % Realizaciones

% Crea las señales
var_r2 = 5e-4;
h_r = [1 0.9 0.5 0.45 0.35 0.25];
h_v = [0.8 0.2 -0.1];

r = normrnd(0, sqrt(var_r2), N, m);
s = filter(h_r,1,r);

sigma_s2 = var(s);
sigma_v2 = sigma_s2 / 100; % SNR = 20db

v = zeros(N,m);
for i=1:m
  v(:,i) = normrnd(0,sqrt(sigma_v2(i)),N,1);  
end

u = filter(h_v, 1, v);
x = s + v;

%% ---------- PARTE C -----------
M = 3;              % Orden del filtro
w0 = ones(M,1)*5;   % Condiciones iniciales
mu = 50;            % Paso

[s_somb,w,error] = canceladorRGB(x,u,s,w0,mu,M);

% Coeficientes del filtro en funcion de las iteraciones
figure(1); hold on;
for i = 1:length(h_v)
  plot(1:1:N, reshape(w(i,1, :),[1,N]))
  line([1,N], [h_v(i), h_v(i)])
end
xlabel('n');
ylabel('$\vec{w}_n$','Interpreter','latex');
title('Coeficientes estimados del filtro');

if m > 1
  J = mean(power(abs(s_somb'),2));
  E = mean(power(abs(error'),2));
else
  J = power(abs(s_somb'),2);
  E = power(abs(error'),2);
end

% Curva de aprendizaje
figure(2); hold on;
plot(1:1:N, J);
line([1 N], [mean(sigma_s2) mean(sigma_s2)],'Color','red');
xlabel('n');
ylabel('$\hat{J}(n)$','Interpreter','latex');
title('Curva de aprendizaje');

% Error
figure(3); hold on;
plot(1:1:N, E);
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');
title('Error de estimacion');

%% ---------- PARTE D -----------

mu = 50;
figure(4); hold on;
% Error en infinito para orden M=1,...,5
for M = 1:5
  w0 = ones(M,1)*5;
  [~,~,error] = canceladorRGB(x,u,s,w0,mu,M);
  E = mean(power(abs(error'),2));
  
  error_inf = mean(E(15000:end));
  stem(M, error_inf);
end
title('Error en infinito vs M');
xlabel("M");
ylabel('$\hat{E}(\infty)$','Interpreter','latex');


%% ---------- PARTE E -----------

m = 50;
M = 2;
w0 = ones(M,1)*5;

figure(5); hold on;
% Error en infinito para pasos mu=10,20,...,100
errores_vs_mu = zeros(8,N);
for mu = 30:10:100

  [~,~,error] = canceladorRGB(x,u,s,w0,mu,M);
  E = mean(power(abs(error'),2));
  errores_vs_mu(mu/10-2,:) = E;
  error_inf = mean(E(15000:end));
  stem(mu, error_inf);
end
title('Error en infinito vs $\mu$','Interpreter','latex');
xlabel('$\mu$','Interpreter','latex');
ylabel('$\hat{E}(\infty)$','Interpreter','latex');

% Error en funcion del paso
figure(6);
semilogy(1:1:N, errores_vs_mu');
title('Error segun $\mu$','Interpreter','latex');
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');
legend('30', '40', '50', '60', '70', '80', '90', '100');

%% ---------- PARTE F -----------

[s, fs] = audioread("Pista_01.wav");

% Ajuste de amplitud
modulacion = sqrt(mean(sigma_s2) / var(s));
s_modulada = modulacion * s;

% Crea ruido y se lo agrega a la señal
N = length(s);
sigma_v2 = mean(sigma_s2) / 100;
v = normrnd(0, sqrt(sigma_v2), N, 1);
u = filter(h_v, 1, v);
x = s + v;

mu = 50;
M = 3;
w0 = ones(M,1)*5;

[s_somb,~,error] = canceladorRGB(x,u,s,w0,mu,M);

E = power(abs(error), 2);

figure(7);
plot(1:length(E), E);
title('Error de estimacion para "Pista_01.wav"');
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');

audiowrite('Pista_01_ruidosa.wav', x, fs);
audiowrite('Pista_01_filtrada.wav', s_somb, fs);


%% PROBLEMA 2

N = 10000;
m = 500;

fs = 44100;
w0 = 2 * pi * 500 / fs;
var_r2 = 5e-4;

A = normrnd(0.1, sqrt(0.003), 1, m);
fase = unifrnd(0,2*pi,1,m);
n = (1:1:N)';

g = A .* sin(w0*n + fase);
h_r = [1 0.9 0.5 0.45 0.35 0.25];
u = [cos(w0*n)' ;  sin(w0*n)'];
r = normrnd(0, sqrt(var_r2), N, m);
s = filter(h_r, 1, r);
x = s + g;

%% -------- PARTE C --------

M = 2;
mu = 1e-3;
w0 = zeros(M,1);

[s_somb,w,error] = canceladorInterferencia(x,u,s,w0,mu,M);

figure(8); hold on;
for i = 1:M
  plot(1:1:N, reshape(w(i,1, :),[1,N]))
end
xlabel('n');
ylabel('$\vec{w}_n$','Interpreter','latex');
title('Coeficientes estimados del filtro');

if m > 1
  J = mean(power(abs(s_somb'),2));
  E = mean(power(abs(error'),2));
else
  J = power(abs(s_somb'),2);
  E = power(abs(error'),2);
end

figure(9); hold on;
plot(1:1:length(J), J);
line([1, N], [mean(var(s)), mean(var(s))]);
xlabel('n');
ylabel('$\hat{J}(n)$','Interpreter','latex');
title('Curva de aprendizaje');

figure(10); hold on;
plot(1:1:length(E), E);
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');
title('Error de estimacion');

%% ------------ PARTE D --------------

M = 2;
w0 = zeros(M,1);

figure(11); hold on;
for mu = 1e-3:1e-3:5e-3

  [~,~,error] = canceladorInterferencia(x,u,s,w0,mu,M);
  E = mean(power(abs(error'),2));

  error_inf = mean(E(end-5000:end));
  stem(mu, error_inf);
end
title('Error en infinito vs $\mu$','Interpreter','latex');
xlabel('$\mu$','Interpreter','latex');
ylabel('$\hat{E}(\infty)$','Interpreter','latex');

%% ------------ PARTE E --------------

N = 10000;
var_r2 = 5e-4;
h_r = [1 0.9 0.5 0.45 0.35 0.25];
r = normrnd(0, sqrt(var_r2), N, 1);
s = filter(h_r, 1, r);
sigma_s2 = var(s);

[s, fs] = audioread("Pista_02.wav");
modulacion = sqrt(mean(sigma_s2) / var(s));
s = modulacion * s;
N = length(s);

A = normrnd(0.1, sqrt(0.003), 1, 1);
fase = unifrnd(0, 2*pi, 1, 1);
n = (1:1:N)';
w0 = 2 * pi * 500 / fs;
g = A * sin(w0*n + fase);
u = [cos(w0*n)' ;  sin(w0*n)'];

x = s + g;

mu = 1e-3;
M = 2;
w0 = zeros(M,1);

[s_somb, ~, error] = canceladorInterferencia(x, u, s, w0, mu,M);
E = power(abs(error), 2);

figure(12);
plot(1:length(E), E);
title('Error de estimacion para "Pista_02.wav');
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');

audiowrite('Pista_02_ruidosa.wav', x, fs);
audiowrite('Pista_02_filtrada.wav', s_somb, fs);

%%  PROBLEMA 3

N = 10000;

var_r2 = 5e-4;
h_r = [1 0.9 0.5 0.45 0.35 0.25];
r = normrnd(0, sqrt(var_r2), N, 1);
s = filter(h_r, 1, r);

sigma_s2 = var(s);
modulacion = sqrt(sigma_s2/var(s));
s_modulada = modulacion * s;

[s, fs] = audioread("Pista_03.wav");

N = length(s);
n = (1:1:N)';

% Crea el ruido blanco
h_v = [0.8 0.2 -0.1];
sigma_v2 = sigma_s2/100;
v = normrnd(0, sqrt(sigma_v2), N, 1);

% Crea interferencia
A = normrnd(0.1,sqrt(0.003), 1, 1);
fase = rand * 2 * pi;
w0 = 2*pi*500/fs;
g = A * sin(w0*n + fase);

% Señales conocidas
x = s + v + g;
u_v = filter(h_v,1,v);
u_g = [cos(w0*n)' ;  sin(w0*n)'];

M = 2;
mu = 2e-3;
w0 = zeros(M,1);
% Aplica LMS para filtrar la interferencia
[sv_somb,~,~] = canceladorInterferencia(x,u_g,s,w0,mu,M);

M = 3;
mu = 10;
w0 = zeros(M,1);
%  Aplica LMS para filtrar el ruido
[~,~,error] = canceladorRGB(sv_somb,u_v,s,w0,mu,M);

% Orden inverso: primero filtrar ruido despues interferencia
% 
% M = 3;
% mu = 10;
% w0 = zeros(M,1);
% %  Aplica LMS para filtrar el ruido
% [sg_somb,~,~] = canceladorRGB(x,u_v,s,w0,mu,M);
% 
% M = 2;
% mu = 2e-3;
% w0 = zeros(M,1);
% % Aplica LMS para filtrar la interferencia
% [~,~,error] = canceladorInterferencia(sg_somb,u_g,s,w0,mu,M);

% error s-ssomb
E = power(abs(error'),2);

figure(13); hold on;
plot(1:1:length(E), E);
title('Error de estimacion para "Pista_03.wav');
xlabel('n');
ylabel('$\hat{E}(n)$','Interpreter','latex');


audiowrite('Pista_03_ruidosa.wav', x, fs);
audiowrite('Pista_03_filtrada.wav', s_somb, fs);

%% FUNCIONES

% Calcula la señal estimada sin ruido balnco
% Los coeficientes del filtro en funcion de las iteraciones y el error, 
% para todas las realizaciones
function [s_somb,w,error] = canceladorRGB(x,u,s,w0,mu,M)
    [N m] = size(x);
    w = zeros(M, m, N);
    w(:,:,1) = ones(M,m).*w0; % Condicion inicial
    v_somb = zeros(N,m);
    s_somb = zeros(N,m);
    error = zeros(N,m);
    for i = M:N
        wi = w(:, :, i-M+1);    % todos los coeficientes de todas las realizaciones en el instante i
        ui = u(i:-1:i-M+1,:);    % una ventana correespondiente al instante i para todas las realizaciones

        v_somb(i,:) = diag(wi' * ui);   % estimacion del ruido para todas las realizaciones en el instante i
        s_somb(i,:) = x(i,:) - v_somb(i,:);

        error(i, :) = s_somb(i,:) - s(i,:);

        w(:, :, i-M+2) = wi + mu * ui .* (ones(M,1)*conj(s_somb(i,:)));

    end
end

% Calcula la señal estimada sin interferencia 
% Los coeficientes del filtro en funcion de las iteraciones y el error, 
% para todas las realizaciones
function [s_somb,w,error] = canceladorInterferencia(x,u,s,w0,mu,M)
    [N m] = size(x);
    w = zeros(M, m, N);
    w(:,:,1) = ones(M,m).*w0; % Condicion inicial
    v_somb = zeros(N,m);
    s_somb = zeros(N,m);
    error = zeros(N,m);
    for i = M:N
        wi = w(:, :, i-M+1);    % todos los coeficientes de todas las realizaciones en el instante i
        ui = u(:,i);            % una ventana correespondiente al instante i para todas las realizaciones
        
        v_somb(i,:) = wi'*ui;   % estimacion del ruido para todas las realizaciones en el instante i
        s_somb(i,:) = x(i,:) - v_somb(i,:);

        error(i, :) = s_somb(i,:) - s(i,:);

        w(:, :, i-M+2) = wi + mu * ui .* (conj(s_somb(i,:)));

    end
end