function z = crtbp (~, y)

% first order circular restricted three-body x-y equations of motion

% input

%  y(1) = normalized x-component of position
%  y(2) = normalized x-component of velocity
%  y(3) = normalized y-component of position
%  y(4) = normalized y-component of velocity

% output

%  z(1) = normalized x component of velocity
%  z(2) = normalized x component of acceleration
%  z(3) = normalized y component of velocity
%  z(4) = normalized y component of acceleration

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu

r1 = sqrt((y(1) + mu) ^ 2 + y(3) ^ 2);

r2 = sqrt((y(1) - 1.0 + mu) ^ 2 + y(3) ^ 2);

z(1) = y(2);

z(2) = 2.0 * y(4) + y(1) - (1.0 - mu) * (y(1) + mu) / r1 ^ 3 - mu ...
    * (y(1) - 1.0 + mu) / r2 ^ 3;

z(3) = y(4);

z(4) = -2.0 * y(2) + y(3) - (1.0 - mu) * y(3) / r1 ^ 3 - mu ...
    * y(3) / r2 ^ 3;
