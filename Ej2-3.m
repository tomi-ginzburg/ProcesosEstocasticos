
clc
clear all
close all

function [y U] = comprimir(x,V,ux,CR,length,width,length_img,width_img)
  datos_totales = length_img*width_img;
  datos_comp = double(CR*datos_totales); % subsize*(length*width+columns(x))+length*width
  subsize = (datos_comp-length*width)/(length*width+columns(x));
  %subsize = double(length*width*CR);
  U = V(:,1:subsize);
  y = U'*(x-ux);
end


img_rgb = imread('img_03.jpg');
img_gray = rgb2gray(img_rgb);
img = im2double(img_gray)*255;
length_img = size(img)(1);
width_img = size(img)(2);

length = 8;
width = 8;
L = length_img*width_img/(length*width);

x = [];
for i=1:fix(width_img/width)
  for j=1:fix(length_img/length)
    x = [x,img((j-1)*length+1:j*length,(i-1)*width+1:i*width)(:)];
  endfor
endfor

ux = mean(x,2);
Cx = cov(x');
[V D] = svd(Cx);

#RECONTRUCCION
[y U] = comprimir(x,V,ux,0.2,length,width,length_img,width_img);
xr = U*y+ux;

img_reconstruida = [];
aux = [];
for i=1:fix(width_img/width)
  for j=1:fix(length_img/length)
    aux = [aux;reshape(xr(:,(i-1)*fix(length_img/length)+j),[length width])];
  endfor
  img_reconstruida = [img_reconstruida aux];
  aux = [];
endfor

img_reconstruida = uint8(img_reconstruida);
figure();
#subplot(1,2,1);
imshow(img_gray);
print('Imagen3 original.png')
figure();
#subplot(1,2,2);
imshow(img_reconstruida);
print('Imagen3 reconstruida.png')

datosTotales = rows(img)*columns(img)
datosComp = rows(y)*columns(y)+rows(ux)+rows(U)*columns(U)
tasaCompresion = datosComp/datosTotales
#tasaCompresion = (rows(y)*columns(y)+rows(ux)+rows(U)*columns(U))/(rows(img)*columns(img))
MSE = sum(((x-xr).*(x-xr))(:))/(rows(x)*columns(x))
