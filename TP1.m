%% EJERCICIO 1
clear all; close all; clc

% Convierte imagenes a escala de grises
img1_rgb = imread('img_01.jpg');
img1_gray = rgb2gray(img1_rgb);
img1 = im2double(img1_gray)*255;

img2_rgb = imread('img_02.jpg');
img2_gray = rgb2gray(img2_rgb);
img2 = im2double(img2_gray)*255;

length_img1 = size(img1,1);
width_img1 = size(img1,2);
length_img2 = size(img2,1);
width_img2 = size(img2,1);

% divide en vectores de 2x1
length = 2;
width = 1;
x1 = dividir_bloques(img1,length,width);
x2 = dividir_bloques(img2,length,width);

corr_x1 = corrcoef(x1');
corr_x2 = corrcoef(x2');

figure(1);
subplot(2,2,1)
imshow(img1_gray);
title('Imagen 1')
print('img_01-Gris.png')

subplot(2,2,2)
imshow(img2_gray);
title('Imagen 2')
print('img_02-Gris.png')

subplot(2,2,3)
scatter(x1(1,:),x1(2,:));
xlabel ("X_0");
ylabel ("X_1");
title('Disperison imagen 1')
print('Dispersion_img1.png')


subplot(2,2,4)
scatter(x2(1,:),x2(2,:));
xlabel ("X_0");
ylabel ("X_1");
title('Disperison imagen 2')
print('Dispersion_img2.png')


%% EJERCICIO 2

img3_rgb = imread('img_03.jpg');
img3_gray = rgb2gray(img3_rgb);
img3 = im2double(img3_gray)*255;

length_img3 = size(img3,1);
width_img3 = size(img3,2);

length = 8;
width = 8;
CR = 0.2;
[y U ux] = comprimir(img3,CR,length,width);

datos_originales = numel(img3);
datos_guardados = numel(y) + numel(U) + numel(ux) ;
CR_real = datos_guardados/datos_originales;

%% EJERCICIO 3

img3_reconstruida = descomprimir(y,U,ux,length,width,length_img3,width_img3);

figure(2);
subplot(1,2,1);
imshow(img3_gray)
title('Original');
subplot(1,2,2);
imshow(img3_reconstruida);
title('Reconstruccion CR=20%')

%% EJERCICIO 4
img4_rgb = imread('img_04.jpg');
img4_gray = rgb2gray(img4_rgb);
img4 = im2double(img4_gray)*255;

length_img4 = size(img4,1);
width_img4 = size(img4,2);

figure(3);
subplot(3,3,1);
imshow(img4_gray);
title("Original");

MSEs = zeros(1,19);
CRs = 5:5:95;
for CR = CRs
  [y U ux] = comprimir(img4,CR/100,length,width);
  img4_reconstruida = descomprimir(y,U,ux,length,width,length_img4,width_img4);
  MSEs(CR/5) = calc_MSE(img4_gray, img4_reconstruida);

  if (CR <= 25)
    k = double(CR/5);
    subplot(3,2,6-k+1);
    imshow(img4_reconstruida);
    title(sprintf("CR = %d %%", CR));
  end
end

figure(4);
stem(CRs, MSEs);
xlabel("Porcentaje de compresión");
ylabel("MSE");


%% FUNCIONES

function x = dividir_bloques(img,length,width)
    length_img = size(img,1);
    width_img = size(img,2);
    L = floor(length_img*width_img/(length*width)); % cantidad de bloques

    x = zeros(width*length,L);
    rows = floor(length_img/length);
    cols = floor(width_img/width);
    for i=1:rows
      for j=1:cols
        bloque = img((i-1)*length+1:i*length,(j-1)*width+1:j*width); 
        x(:, j + (i-1)*cols) = bloque(:);
      end
    end
end

function [y U ux] = comprimir(img,CR,length,width)

    length_img = size(img,1);
    width_img = size(img,2);
    
    x = dividir_bloques(img,length,width);
    ux = mean(x,2);
    Cx = cov(x');
    
    [V D] = svd(Cx);
    
    datos_totales = length_img*width_img;
    datos_comp = double(CR*datos_totales);
   
    % Subsize: Cálculo de cuantas columnas quedarse de la matriz de autovalores
    % Si length=8 y width=8 -> bloques de 8x8
    % U = 64 x subsize , Y = subsize x columnas(X) 
    subsize = floor((datos_comp-size(ux,1))/(length*width+size(x,2)));
    U = V(:,1:subsize);
    y = U'*(x-ux);
end

function img_reconstruida = descomprimir(y,U,ux,length,width,length_img,width_img)
    xr = U*y+ux;
    cols = floor(width_img/width);
    rows = floor(length_img/length);
    img_reconstruida = zeros(rows*length,cols*width);
    for i=1:rows
      for j=1:cols
        bloque = reshape(xr(:,(i-1)*cols+j),[length width]);
        img_reconstruida(length*(i-1) + 1 : length*i, width*(j-1) + 1 : width*j) = bloque;
      end
    end

    img_reconstruida = uint8(img_reconstruida);
end

function MSE = calc_MSE(image, reconstruido)
    length_img = size(reconstruido,1);
    width_img = size(reconstruido,2);
    error = image(1:length_img,1:width_img)-reconstruido;
    MSE = sum(error.^2,'all')/(length_img*width_img);
end