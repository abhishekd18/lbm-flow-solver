% Load the data files
load Ux.dat
load Uy.dat
load ny.dat
load nx.dat

% Geometry representation
Nx = 70;
Ny = 70;
rsqr = 0.05;
rx = sqrt(rsqr)*Nx;
ry = sqrt(rsqr)*Ny;
theta = [0:pi/100:2*pi];

X = rx*cos(theta) + Nx/2;
Y = ry*sin(theta) + Ny/2;

% Plot contours
f=figure(1);
set(f, 'Position', [0 0 850 800]);
plot(X,Y,'black','linewidth',2);
axis([1 Nx 1 Ny]);
axis equal;
title('Velocity Magnitude |U| contours');
xlabel('Nx');
ylabel('Ny');
hold on;
contourf(nx, ny, sqrt(Ux.^2+Uy.^2),'ShowText','on');
saveas(f,'contour.png','png');
close(f);

% Plot vectors
f=figure(1);
set(f, 'Position', [0 0 1000 1000]);
quiverc(nx, ny, Ux, Uy, 2);
hold on;
plot(X,Y,'black','linewidth',2);
axis([1 Nx 1 Ny]);
axis equal;
title('Velocity Vectors');
xlabel('Nx');
ylabel('Ny');
saveas(f,'vectors.png','png');
close(f);

% Streamlines
streamline(Ux,Uy,startx',starty);
