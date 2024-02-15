clear all; close all; clc

%% EJERCICIO 2

files2 = char('audio_a.wav', 'audio_e.wav', 'audio_s.wav', 'audio_sh.wav');

% Orden
p = [5, 10, 30];

% Window time en seg (30ms)
Tw = 30e-3;
nfft = 1024;

for i = 1:size(files2,1)
   % Muestras, frecuencia de muestreo en Hz
  [y Fs] = audioread(files2(i, :));
  len_y = length(y);

  % Ventana hamming
  cantidad_muestras  = round(Fs*Tw);
  window = hamming(cantidad_muestras);

  % Obtener una ventana en el medio de la señal
  xs = y(round((len_y - cantidad_muestras)/2):round((len_y + cantidad_muestras)/2)-1);
  % Suavizado
  xs = xs .* window;

  figure(i);

  for j = 1:length(p)

    % Parámetros LPC de la ventana elegida
    [a G] = param_lpc(xs,p(j));

    % PSD estimada = abs(H).^2 * S_x
    w = linspace(0, 2*pi, nfft);
    S_estimada = psd(G, a, w);

    % Periodograma sesgado
    periodograma = (1 / length(xs)) * abs(fft(xs, nfft)) .^ 2;

    % Autocorrelación
    autocor = xcorr(xs);

    % Gráficos
    % Respuesta temporal
    s1 = subplot(3, 1, 1);
    plot(xs);
    ylabel(sprintf ("Respuesta temporal"), 'FontWeight','bold');

    % Autocorrelacion
    s2 = subplot(3, 1, 2);
    len = (length(autocor) - 1) / 2;
    plot(-len:1:len, autocor);
    ylabel (sprintf ("Autocorrelacion"), 'FontWeight','bold');

    % PSD
    s3 = subplot(3, 3, j + 6); hold on;
    plot(w, 20*log10(S_estimada));
    plot(w, 20*log10(periodograma));
    legend('Analítico','Periodograma');
    if j==1
      ylabel (sprintf ("PSD"), 'FontWeight','bold');
    end
    xlabel (sprintf ("p = %d", p(j)), 'FontWeight','bold');
    legend ("location", "southeast")

  end
  sgtitle(files2(i,:));
end

%% EJERCICIO 3

% . . . . . . . . . . . .
% ------- ------- -------
%     ------- -------

% Window time en seg (30ms)
Tw = 50e-3;
% Orden
p = 15;
% Umbral
alpha = 0.15;
% Compensación de amplitud para la reconstrucción
compensacion = 5;

for i = 1:size(files2,1)

  [y Fs] = audioread(files2(i,:));
  len_y = length(y);

  cantidad_muestras  = round(Fs*Tw);

  % Ventana
  window = hamming(cantidad_muestras);
  % Pitch fijo para la reconstrucción [Hz]
  if (i<=2)
    pitch_fijo = 200;
  else
      pitch_fijo = 0;
  end
  % Obtención de parámetros para cada segmento
  [cantidad_segmentos, mat_a, vec_G] = param_lpc_full(y, cantidad_muestras, window, p);

  % Reconstrucción
  x_reconstruida = reconstruir(cantidad_segmentos, cantidad_muestras, alpha, window, y, mat_a, vec_G, Fs, compensacion, pitch_fijo);

  % sound([y', zeros(1, 1000), x_reconstruida], Fs);
  if i==4
    filename = sprintf ("%s_lpc.wav",files2(i, 1:end-4));
  else
    filename = sprintf ("%s_lpc.wav",files2(i, 1:end-5));
  end
  audiowrite(filename, x_reconstruida, Fs)
end

%% EJERCICIO 4 

% Parámetros:
% usa los mismos que el ejercicio 3
files4 = char('audio_01.wav', 'audio_02.wav', 'audio_03.wav', 'audio_04.wav');

for i = 1:size(files4,1)

  [y Fs] = audioread(files4(i, :));

  % Ventana hamming
  cantidad_muestras  = round(Fs*Tw);
  window = hamming(cantidad_muestras);

  [cantidad_segmentos, mat_a, vec_G] = param_lpc_full(y, cantidad_muestras, window, p);

  % Reconstrucción
  x_reconstruida = reconstruir(cantidad_segmentos, cantidad_muestras, alpha, window, y, mat_a, vec_G, Fs, compensacion, -1);

  % sound([y', zeros(1, 1000), x_reconstruida], Fs);
  filename = sprintf ("%s_lpc.wav",files4(i, 1:end-4));
  audiowrite (filename, x_reconstruida, Fs)

end

%%  ANALISIS EXTRA 

% Análisis para la muestra del ejercicio 1
e = filter([1 -a'], 1, xs);
re = xcorr(e)/length(e);
re_truncado = re((length(re) + 1)/2 + 1:end);
[max_secundario, max_secundario_idx] = max(re_truncado);
f0 = Fs/max_secundario_idx;

figure(); hold on;
plot(1:length(re), re / max(re));
line([1, length(re)], [0.1, 0.1], 'Color', 'r');
line([(length(re) + 1)/2 + max_secundario_idx, (length(re) + 1)/2 + 1 + max_secundario_idx], [0, 1], 'Color', 'r');
title("Cálculo de pitch (autocorrelación del residuo)");

% Análisis para un segmento intermedio
figure(); hold on
nfft = 1024;
plot(1:nfft, abs(fft(y(4*cantidad_muestras:5*cantidad_muestras),nfft)));
plot(1:nfft, abs(fft(x_reconstruida(4*cantidad_muestras:5*cantidad_muestras),nfft)));
legend('FFT original', 'reconstrucción');

% Analisis temporal de los audios
figure();
s = subplot(3, 1, 1);
plot(y);
title('Audio original') 
s = subplot(3, 1, 2);
plot(x_reconstruida);
title('Audio reconstruido') 
s = subplot(3, 1, 3);
plot(y(1:length(x_reconstruida'))-x_reconstruida');
title('Diferencia') 

%% Funciones

% Estima los parametros LPC(G,ak) de xs con orden P
function [a G] = param_lpc(xs, P)
  % vector de correlación, orden P (mide 2P + 1)
  r = xcorr(xs, P) / length(xs);

  R = zeros(P,P);
  % agrega fila por fila
  for i = 0:P-1
    % para i = 0, start_idx debe ser r(0), índice P+1
    start_idx = P + 1 - i;
    % para i = 0, end_idx debe ser r(p-1), índice 2P
    end_idx = 2*P - i;
    R(i+1,:) = r(start_idx:end_idx);
  end

  auto_corr = r(P+1); % Valor de r en 0
  r_corto = r(P+2:2*P+1); % Va de 1 hasta P

  a = inv(R) * r_corto;
  G = sqrt(auto_corr - a' * r_corto);

end

%Recibe w como vector de frecuencias
function [mod_H_cuad] = psd(G, a, w)
  sumatoria = zeros(1, length(w));
  % p = length(a)
  for k = 1:length(a)
    sumatoria = sumatoria + a(k) * exp(-1i*k*w);
  end
  mod_H_cuad = abs(G^2 * power(1 - sumatoria, -2));
end

% Genera u de cantidad_muestras de largo
% Si frec_reconstruccion es 0 con ruido haussiano blanco
% si es distinto de 0 genera un tren de impulsos con esa frecuencia
function u = generar_u(frec_reconstruccion, cantidad_muestras, Fs)

  if frec_reconstruccion == 0
    % Usamos 0.4 para ajustar amplitud y que no tenga más volumen que los trenes de deltas
    u = normrnd(0, 0.4, 1, cantidad_muestras);
    return;
  end

  % Vocales
  % Un periodo a la frecuencia de reconstrucción
  cantidad_ceros = round(Fs/frec_reconstruccion) - 1;
  periodo_u = [1 zeros(1, cantidad_ceros)];
  % Cantidad de veces que necesitamos ese periodo y cantidad de muestras sobrantes
  reps_u = floor(cantidad_muestras/(cantidad_ceros+1));
  resto_u = mod(cantidad_muestras, (cantidad_ceros+1));
  % Señal de entrada para reconstrucción de vocales
  % (Tren de impulsos a la frecuencia dada)
  u = repmat(periodo_u, 1, reps_u);
  u = [u periodo_u(1:resto_u)];
end

% Estima el pitch de la señal xs a partir de la autocorrelacion del error
% Si el segundo pico de correlacion no supera el umbral es ruido blanco y
% devuelve f0 = 0.
function f0 = pitch_lpc(xs, a, alpha, Fs)

  % Residuo (del proceso autorregresivo de la salida del filtro)
  e = filter([1 -a'], 1, xs);
  % Autocorrelación del residuo
  re = xcorr(e)/length(e);
  re = re / max(re);
  
  % Mitad derecha de la autocorrelación (sin contar 0)
  re_truncado = re((length(re) + 1)/2 + 1:end);

  % Pico del pitch
  [max_secundario max_secundario_idx] = max(re_truncado);

  % Verificar umbral
  if max_secundario  < alpha
    % Simboliza ruido blanco
    f0 = 0;
    return;
  end

  % Fórmula del PDF
  f0 = Fs / max_secundario_idx;

end

% Recontruye la señal completa. Genera la entrada u dependiendo si usa un
% pitch fijo para todo el audio o lo calcula para cada ventana, y utiliza 
% los coeficientes del filtro para estimar cada ventana de señal y las une
function x_reconstruida = reconstruir(cantidad_segmentos, cantidad_muestras, alpha, window, y, mat_a, vec_G, Fs, compensacion, pitch_fijo)

  cant_total_muestras = ceil(cantidad_segmentos/2)*cantidad_muestras + mod(cantidad_segmentos+1, 2) * floor(cantidad_muestras/2);

  x_reconstruida = zeros(1, cant_total_muestras);
  for i = 1:cantidad_segmentos
    start_idx = (i-1) * floor(cantidad_muestras/2) + 1;
    end_idx = start_idx + cantidad_muestras - 1;

    xs = y(start_idx:end_idx) .* window;

    if pitch_fijo == -1
      f0 = pitch_lpc(xs, mat_a(:, i), alpha, Fs);
    else
      f0 = pitch_fijo;
    end

    u = generar_u(f0, cantidad_muestras, Fs);
    xr = filter(vec_G(i), [1, - mat_a(:,i)'], u);
    x_reconstruida(start_idx:end_idx) = x_reconstruida(start_idx:end_idx) + xr;

  end

  %Ajuste de amplitud
  x_reconstruida = x_reconstruida * compensacion;

end

%Obtiene los parámetros para TODO el audio dado, utilizando
%la cantidad de muestras por ventana dada, con el orden P
function [cantidad_segmentos, mat_a, vec_G] = param_lpc_full(y, cantidad_muestras, window, p)

  cantidad_segmentos = floor(length(y)*2/cantidad_muestras) - 1;
  mat_a = [];
  vec_G = [];

  for i=1:cantidad_segmentos
    start_idx = (i-1) * ceil(cantidad_muestras/2) + 1;
    end_idx = start_idx + cantidad_muestras - 1;
    xs = y(start_idx:end_idx) .* window;
    [a G] = param_lpc(xs, p);
    mat_a = [mat_a a]; %Cada instante/ventana está en cada columna
    vec_G = [vec_G G];
  end

end
