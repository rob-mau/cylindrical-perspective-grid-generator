% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Generatore Griglia per Prospettiva Cilindrica (Disegno)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all
close all

% parametri griglia cubica
% la griglia cubica è costituita da 3 famiglie di segmenti tra loro
% paralleli ed equidistanziati (distanza di 1 "unità"). Le famiglie sono
% tra loro mutuamente ortogonali. Ogni segmento è rappresentato con
% "num_points" punti.
half_len = 2;       % semi-lunghezza del lato del cubo ("unità" intere!)
num_points = 400;   % numero di punti in ciascun segmento

% parametri geometrici
% posizione del centro della griglia cubica
d = [0; 0; 0];  
% rotazione che descrive l'orientamento della griglia (param. asse-angolo)
n = [1; 1; 1];  % asse della rotazione 
angle = pi/8; 	% angolo della rotazione
% raggio della superficie cilindrica
radius = 1; 

% matrice di trasformazione T tra due sistemi di riferimento genericamente
% orientati (matrice omogenea)
n = n / norm(n);
I = eye(3);
A = n * n.';    % proiettore ortogonale sullo span{n}
N = skew(n);   % matrice antisimmetrica generata da n 
% formula di Rodrigues
R = A + (I - A) * cos(angle) + N * sin(angle);
% matrice omogenea che codifica la roto-traslazione
T = [R, d; [0, 0, 0, 1]];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Inizializzazione Matrici
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% le lettere x, y, z individuano ciascuna delle 3 famiglie di 
% segmenti mutuamente ortogonali che compongono la griglia cubica
% (allineate, rispettivamente, lungo gli assi x, y e z)

max_index  = 2 * half_len + 1;

% famiglie di segmenti allineate con gli assi coordinati (pedice a)
% inizializzazione
px_a = zeros(max_index, max_index, 4, num_points);
py_a = zeros(max_index, max_index, 4, num_points);
pz_a = zeros(max_index, max_index, 4, num_points);

% famiglie di segmenti roto-traslate tramite la trasform. T (pedice b)
% inizializzazione
px_b = zeros(max_index, max_index, 4, num_points);
py_b = zeros(max_index, max_index, 4, num_points);
pz_b = zeros(max_index, max_index, 4, num_points);

% famiglie di segmenti proiettate sul cilindro --> diventano tre famiglie
% di curve, in particolare archi di ellissi (pedice c)
% inizializzazione
px_c = zeros(max_index, max_index, 3, num_points);
py_c = zeros(max_index, max_index, 3, num_points);
pz_c = zeros(max_index, max_index, 3, num_points);

% famiglie di curve sul cilindro "srotolato" --> gli archi di ellissi
% diventano archi di sinusoidi (pedice d)
% inizializzazione
px_d = zeros(max_index, max_index, 2, num_points);
py_d = zeros(max_index, max_index, 2, num_points);
pz_d = zeros(max_index, max_index, 2, num_points);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Calcolo Matrici
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% famiglie di segmenti allineate con gli assi coordinati (pedice a)
step = 2 * half_len / (num_points - 1);

x = -half_len:step:half_len;
y = -half_len:step:half_len;
z = -half_len:step:half_len;

for i = 1:max_index
    for j = 1:max_index
        
        px_a(i,j,:,:) = [x; 
                         ones(1, num_points) * ((i - 1) - half_len); 
                         ones(1, num_points) * ((j - 1) - half_len); 
                         ones(1, num_points)];

        py_a(i,j,:,:) = [ones(1, num_points) * ((j - 1) - half_len);
                         y;
                         ones(1, num_points) * ((i - 1) - half_len); 
                         ones(1, num_points)];

        pz_a(i,j,:,:) = [ones(1, num_points) * ((i - 1) - half_len); 
                         ones(1, num_points) * ((j - 1) - half_len);
                         z;
                         ones(1, num_points)];
    end
end
     
% famiglie di segmenti roto-traslate tramite la trasform. T (pedice b)
for i = 1:max_index
    for j = 1:max_index
        
        px_b(i,j,:,:) = T * squeeze(px_a(i,j,:,:));
        py_b(i,j,:,:) = T * squeeze(py_a(i,j,:,:));
        pz_b(i,j,:,:) = T * squeeze(pz_a(i,j,:,:));
    end
end

% famiglie di segmenti proiettate sul cilindro --> diventano tre famiglie
% di curve, in particolare archi di ellissi (pedice c)
for i = 1:max_index
    for j = 1:max_index
        
        px_c(i,j,:,:) = radius/sqrt(px_b(i,j,1,:).^2 + ...
                        px_b(i,j,2,:).^2) .* px_b(i,j,1:3,:);
        py_c(i,j,:,:) = radius/sqrt(py_b(i,j,1,:).^2 + ...
                        py_b(i,j,2,:).^2) .* py_b(i,j,1:3,:);
        pz_c(i,j,:,:) = radius/sqrt(pz_b(i,j,1,:).^2 + ...
                        pz_b(i,j,2,:).^2) .* pz_b(i,j,1:3,:);
    end
end

% famiglie di curve sul cilindro "srotolato" --> gli archi di ellissi
% diventano archi di sinusoidi (pedice d)
for i = 1:max_index
    for j = 1:max_index
        
        px_d(i,j,1,:) =  radius * atan2(px_c(i,j,2,:), px_c(i,j,1,:));
        px_d(i,j,2,:) =  px_c(i,j,3,:);
        
        py_d(i,j,1,:) =  radius * atan2(py_c(i,j,2,:), py_c(i,j,1,:));
        py_d(i,j,2,:) =  py_c(i,j,3,:);
        
        pz_d(i,j,1,:) =  radius * atan2(pz_c(i,j,2,:), pz_c(i,j,1,:));
        pz_d(i,j,2,:) =  pz_c(i,j,3,:);   
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Figura 1 - Griglia Cubica Tridimensionale & Cilindro
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f1 = figure('Name', 'Cilindro', 'Color', 'w');
axis vis3d equal off
axis([-1, 1, -1, 1, -1, 1] * 10);
hold on

% superficie cilindrica
[t, h] = meshgrid(0:pi/10:2*pi, -8:0.5:8);
xc = radius .* cos(t);
yc = radius .* sin(t);
zc = h;
surf(xc, yc, zc, 'FaceColor', 'k', 'EdgeColor', 'none')
alpha 0.05

% centro del sistema di riferimento e centro della griglia cubica
plot3(0, 0, 0, '*')
plot3(d(1), d(2), d(3), '*')

for i = 1:max_index
    for j = 1:max_index
        % segmenti griglia cubica
        plot3(squeeze(px_b(i,j,1,:)), ...
              squeeze(px_b(i,j,2,:)), ...
              squeeze(px_b(i,j,3,:)), 'Color', 'b')
        plot3(squeeze(py_b(i,j,1,:)), ...
              squeeze(py_b(i,j,2,:)), ...
              squeeze(py_b(i,j,3,:)), 'Color', 'c')
        plot3(squeeze(pz_b(i,j,1,:)), ...
              squeeze(pz_b(i,j,2,:)), ...
              squeeze(pz_b(i,j,3,:)), 'Color', 'm')
        
        % proiezioni sulla superf. cilindrica (archi di ellisse)
        plot3(squeeze(px_c(i,j,1,:)), ...
              squeeze(px_c(i,j,2,:)), ...
              squeeze(px_c(i,j,3,:)), 'Color', 'b')
        plot3(squeeze(py_c(i,j,1,:)), ...
              squeeze(py_c(i,j,2,:)), ...
              squeeze(py_c(i,j,3,:)), 'Color', 'c')
        plot3(squeeze(pz_c(i,j,1,:)), ...
              squeeze(pz_c(i,j,2,:)), ... 
              squeeze(pz_c(i,j,3,:)), 'Color', 'm')
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Figura 2 - Griglia Prospettica Bidimensionale
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f2 = figure('Name', 'Proiezione', 'Color', 'w');
axis equal
axis([-pi, pi, -pi/2, pi/2]);
hold on

% tracciamento curve della griglia prospettica 
for i = 1:max_index
    for j = 1:max_index
        
        plot(-squeeze(pz_d(i,j,1,:)), squeeze(pz_d(i,j,2,:)), ...
             '.', 'MarkerSize', 1, 'Color', 'b')
        plot(-squeeze(py_d(i,j,1,:)), squeeze(py_d(i,j,2,:)), ...
             '.', 'MarkerSize', 1, 'Color', 'c')
        plot(-squeeze(px_d(i,j,1,:)), squeeze(px_d(i,j,2,:)), ...
             '.', 'MarkerSize', 1, 'Color', 'm')
    end
end

% tracciamento dei riferimenti
plot(-squeeze(px_d(1,1,1,:)), squeeze(px_d(1,1,2,:)), ...
     '.', 'MarkerSize', 3, 'Color', 'k')
plot(-squeeze(px_d(max_index,max_index,1,:)), ...
      squeeze(px_d(max_index,max_index,2,:)), ... 
     '.', 'MarkerSize', 3, 'Color', 'k')

plot(-squeeze(py_d(1,1,1,:)), squeeze(py_d(1,1,2,:)), ...
     '.', 'MarkerSize', 3, 'Color', 'k')
plot(-squeeze(py_d(max_index,max_index,1,:)), ...
      squeeze(py_d(max_index,max_index,2,:)), ...
     '.', 'MarkerSize', 3, 'Color', 'k')

plot(-squeeze(pz_d(1,1,1,:)), squeeze(pz_d(1,1,2,:)), ...
     '.', 'MarkerSize', 3, 'Color', 'k')
plot(-squeeze(pz_d(max_index,max_index,1,:)), ...
      squeeze(pz_d(max_index,max_index,2,:)), ...
     '.', 'MarkerSize', 3, 'Color', 'k')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Definizione Funzioni
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Funzione out = skew(in):
%   in:  vettore tridimensionale
%   out: matrice antisimmetrica 3*3
function out = skew(in)
    a = in(1);
    b = in(2);
    c = in(3);

    out = [0, -c,  b;
           c,  0, -a;
          -b,  a,  0];
end
