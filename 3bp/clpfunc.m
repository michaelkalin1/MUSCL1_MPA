function y = clpfunc(x)

% collinear libration points
% normalized x-coordinate objective function

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ilp mu

switch ilp
    
    case 1
        
        y = x - (1.0 - mu) / (mu + x)^2 + mu / (x - 1.0 + mu)^2;
        
    case 2
        
        y = x - (1.0 - mu) / (mu + x)^2 - mu / (x - 1.0 + mu)^2;
        
    case 3
        
        y = x + (1.0 - mu) / (mu + x)^2 + mu / (x - 1.0 + mu)^2;
        
end
