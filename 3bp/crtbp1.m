% crtbp1.m        February 11, 2022

% compute x and y coordinates and energy of the
% equilibrium libration points of the circular-restricted
% three-body problem (crtbp). Display geometry graphics.

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global ilp mu

clc; home;

fprintf('\n            program crtbp1\n');

fprintf('\n< equilibrium coordinates and energy >\n\n');

while(1)
    fprintf('\nplease input the value for the mass ratio\n');

    mu = input('? ');

    if (mu > 0)
        break;
    end
end

% location of the larger mass (normalized)

xm1 = -mu;

% location of the smaller mass (normalized)

xm2 = 1.0 - mu;

% L1 libration point (normalized)

ilp = 1;

xr1 = -2;

xr2 = +2;

rtol = 1.0e-8;

[xl1, ~] = brent ('clpfunc', xr1, xr2, rtol);

yl1 = 0;

r1sqr = (xl1 - xm1)^2 + yl1^2;

r2sqr = (xl1 - xm2)^2 + yl1^2;

e1 = -0.5 * (xl1^2 + yl1^2) - (1.0 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

% L2 libration point (normalized)

ilp = 2;

xr1 = -2.0;

xr2 = +2.0;

rtol = 1.0e-8;

[xl2, ~] = brent ('clpfunc', xr1, xr2, rtol);

yl2 = 0.0;

r1sqr = (xl2 - xm1)^2 + yl2^2;

r2sqr = (xl2 - xm2)^2 + yl2^2;

e2 = -0.5 * (xl2^2 + yl2^2) - (1.0 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

% L3 libration point (normalized)

ilp = 3;

xr1 = -2.0;

xr2 = +2.0;

rtol = 1.0e-8;

[xl3, froot] = brent ('clpfunc', xr1, xr2, rtol);

yl3 = 0.0;

r1sqr = (xl3 - xm1)^2 + yl3^2;

r2sqr = (xl3 - xm2)^2 + yl3^2;

e3 = -0.5 * (xl3^2 + yl3^2) - (1.0 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

% L4 libration point (normalized)

xl4 = 0.5 - mu;

yl4 = 0.5 * sqrt(3);

r1sqr = (xl4 - xm1)^2 + yl4^2;

r2sqr = (xl4 - xm2)^2 + yl4^2;

e4 = -0.5 * (xl4^2 + yl4^2) - (1.0 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

% L5 libration point (normalized)

xl5 = 0.5 - mu;

yl5 = - 0.5 * sqrt(3);

r1sqr = (xl5 - xm1)^2 + yl5^2;

r2sqr = (xl5 - xm2)^2 + yl5^2;

e5 = -0.5 * (xl5^2 + yl5^2) - (1.0 - mu) / sqrt(r1sqr) - mu / sqrt(r2sqr);

% print numerical results

fprintf('\n                          program crtbp1\n');

fprintf('\n              < equilibrium coordinates and energy >\n\n');

fprintf('\nmass ratio = %12.10e\n', mu);

fprintf('\nlocation      x-coordinate       y-coordinate            energy\n');

fprintf('\n   L1       %14.10f     %14.10f      %12.10e\n', xl1, yl1, e1);

fprintf('\n   L2       %14.10f     %14.10f      %12.10e\n', xl2, yl2, e2);

fprintf('\n   L3       %14.10f     %14.10f      %12.10e\n', xl3, yl3, e3);

fprintf('\n   L4       %14.10f     %14.10f      %12.10e\n', xl4, yl4, e4);

fprintf('\n   L5       %14.10f     %14.10f      %12.10e\n\n', xl5, yl5, e5);

%%%%%%%%%%%%%%%%%%%%%%%%
% display CR3BP geometry
%%%%%%%%%%%%%%%%%%%%%%%%

mu_str = num2str(mu, 'mu = %12.10e\n');

% default values for plot boundaries

xmin = -1.5;
xmax = +1.5;
ymin = -1.5;
ymax = +1.5;

axis ([xmin xmax ymin ymax]);

axis square;

ylabel('normalized y coordinate', 'FontSize', 14);

xlabel('normalized x coordinate', 'FontSize', 14);

% label locations of larger mass (green) and smaller mass (blue)

hold on;

plot(-mu, 0, 'og', 'MarkerSize', 16, 'MarkerFaceColor', 'g');

plot(1.0 - mu, 0, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% plot and label libration points

plot(xl1, yl1, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

text(xl1 - 0.05, yl1 + 0.1, 'L_1');
    
plot(xl2, yl2, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

text(xl2 - 0.05, yl2 + 0.1, 'L_2');

plot(xl3, yl3, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

text(xl3 - 0.05, yl3 + 0.1, 'L_3');

plot(xl4, yl4, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

text(xl4 - 0.05, yl4 + 0.1, 'L_4');

plot(xl5, yl5, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

text(xl5 - 0.05, yl5 - 0.1, 'L_5');

% display x and y axes

plot([xmin xmax], [0.0 0.0], '-k');

plot([0.0 0.0], [ymin ymax], '-k');

% display equilateral graphics

plot([-mu xl4], [0.0 yl4], '-k');

plot([-mu xl5], [0.0 yl5], '-k');

plot([1.0 - mu xl4], [0.0 yl4], '-k');

plot([1.0 - mu xl5], [0.0 yl5], '-k');

text(-1.4, +1.3, mu_str, 'FontSize', 12);

text(-1.4, +1.1, 'barycenter at (0,0)', 'FontSize', 12);

text(-1.4, +0.9, 'mass sizes not to scale', 'FontSize', 12);

text(-1.4, +0.7, 'major body - green', 'FontSize', 12, 'Color', 'green');

text(-1.4, +0.5, 'minor body - blue', 'FontSize', 12, 'Color', 'blue');

title({'Geometry of the CR3BP', '(rotating coordinate system)'},'FontSize', 16);

% create tiff graphics disk file

print ('cr3bp.tif', '-dtiff'); 

fprintf('\n\n');
