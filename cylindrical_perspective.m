% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Grid Generator for Cylindrical Perspective (Drawing)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; 
close all; 
clc;

% Cubic Grid Parameters
% The cubic grid consists of 3 sets of mutually parallel and equidistant 
% segments (spaced by 1 "unit"). The sets are mutually orthogonal. 
% Each segment is represented by "num_points" points.
half_length = 2;    % half-length of the cube side (integer number)
num_points = 400;   % number of points in each segment

% Geometric Parameters
% position of the center of the cubic grid
position = [0; 0; 0];  
% orientation of the grid (axis-angle parametrization)
rot_axis = [1; 1; 1];  % rotation axis
rot_angle = pi/8; 	   % rotation angle
% radius of the cylindrical surface
radius = 1; 

% transformation matrix T between two arbitrarily oriented frames
% (homogeneous matrix)
S = skew(rot_axis / norm(rot_axis));
I = eye(3);
% Rodrigues' rotation formula
orientation = I + sin(rot_angle) * S + (1 - cos(rot_angle)) * S^2;

% homogeneous matrix encoding the rotation and translation
T = [orientation, position; [0, 0, 0, 1]];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Matrix Initialization
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The letters x, y, z identify each of the 3 families of segments 
% that are mutually orthogonal composing the cubic grid 
% (aligned along the x, y, and z axes, respectively)
max_index  = 2 * half_length + 1;

% families of segments aligned with the coordinate axes (subscript a)
% initialization
px_a = zeros(max_index, max_index, 4, num_points);
py_a = zeros(max_index, max_index, 4, num_points);
pz_a = zeros(max_index, max_index, 4, num_points);

% families of segments transformed through T transformation (subscript b)
% initialization
px_b = zeros(max_index, max_index, 4, num_points);
py_b = zeros(max_index, max_index, 4, num_points);
pz_b = zeros(max_index, max_index, 4, num_points);

% families of segments projected onto the cylinder -- > they become three 
% families of curves, specifically arcs of ellipses. (subscript c)
% initialization
px_c = zeros(max_index, max_index, 3, num_points);
py_c = zeros(max_index, max_index, 3, num_points);
pz_c = zeros(max_index, max_index, 3, num_points);

% families of curves on the 'unrolled' cylinder --> the arcs of ellipses 
% become arcs of sinusoids. (subscript d)
% initialization
px_d = zeros(max_index, max_index, 2, num_points);
py_d = zeros(max_index, max_index, 2, num_points);
pz_d = zeros(max_index, max_index, 2, num_points);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Matrix Computation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% families of segments aligned with the coordinate axes (subscript a)
step = 2 * half_length / (num_points - 1);

x = -half_length:step:half_length;
y = -half_length:step:half_length;
z = -half_length:step:half_length;

for i = 1:max_index
    for j = 1:max_index
        px_a(i,j,:,:) = [x; 
                         ones(1, num_points) * ((i - 1) - half_length); 
                         ones(1, num_points) * ((j - 1) - half_length); 
                         ones(1, num_points)];

        py_a(i,j,:,:) = [ones(1, num_points) * ((j - 1) - half_length);
                         y;
                         ones(1, num_points) * ((i - 1) - half_length); 
                         ones(1, num_points)];

        pz_a(i,j,:,:) = [ones(1, num_points) * ((i - 1) - half_length); 
                         ones(1, num_points) * ((j - 1) - half_length);
                         z;
                         ones(1, num_points)];
    end
end
     
% families of segments transformed through T transformation (subscript b)
for i = 1:max_index
    for j = 1:max_index
        px_b(i,j,:,:) = T * squeeze(px_a(i,j,:,:));
        py_b(i,j,:,:) = T * squeeze(py_a(i,j,:,:));
        pz_b(i,j,:,:) = T * squeeze(pz_a(i,j,:,:));
    end
end

% families of segments projected onto the cylinder -- > they become three 
% families of curves, specifically arcs of ellipses. (subscript c)
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

% families of curves on the 'unrolled' cylinder --> the arcs of ellipses 
% become arcs of sinusoids. (subscript d)
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
%  Figure 1 - 3D Cubic Grid & Cylinder
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_1 = figure('Name', '3D Cubic Grid & Cylinder Projection', 'Color', 'w');
axis vis3d equal off
axis([-1, 1, -1, 1, -1, 1] * 10);
hold on

% cylindrical surface
[t, h] = meshgrid(0:pi/10:2*pi, -8:0.5:8);
xc = radius .* cos(t);
yc = radius .* sin(t);
zc = h;
surf(xc, yc, zc, 'FaceColor', 'k', 'EdgeColor', 'none')
alpha 0.05

% center of the reference system and center of the cubic grid
plot3(0, 0, 0, '*')
plot3(position(1), position(2), position(3), '*')

for i = 1:max_index
    for j = 1:max_index
        % cubic grid segments
        plot3(squeeze(px_b(i,j,1,:)), ...
              squeeze(px_b(i,j,2,:)), ...
              squeeze(px_b(i,j,3,:)), 'Color', 'b')
        plot3(squeeze(py_b(i,j,1,:)), ...
              squeeze(py_b(i,j,2,:)), ...
              squeeze(py_b(i,j,3,:)), 'Color', 'c')
        plot3(squeeze(pz_b(i,j,1,:)), ...
              squeeze(pz_b(i,j,2,:)), ...
              squeeze(pz_b(i,j,3,:)), 'Color', 'm')
        
        % projections onto the cylindrical surface (arcs of ellipses)
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
%  Figura 2 - 2D Perspective Grid
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig_2 = figure('Name', '2D Projection', 'Color', 'w');
axis equal
axis([-pi, pi, -pi/2, pi/2]);
hold on

% drawing curves of the perspective grid
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

% drawing reference frames
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
%  Function Definitions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Funzion out = skew(in):
%   in:  3D vector
%   out: 3*3 skew-symmetric matrix
function out = skew(in)
    a = in(1);
    b = in(2);
    c = in(3);

    out = [0, -c,  b;
           c,  0, -a;
          -b,  a,  0];
end
