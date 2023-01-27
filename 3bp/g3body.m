% g3body.m      February 6, 2022

% graphics display of integrated orbital motion in the
% circular-restricted Earth-Moon three body problem

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu

clc; home;

fprintf('\n             program g3body\n');

fprintf('\n< graphics display of three-body motion >\n\n');

while(1)
    
    fprintf('\n <1> periodic orbit about L1\n\n');

    fprintf(' <2> periodic orbit about L2\n\n');

    fprintf(' <3> periodic orbit about L3\n\n');

    fprintf(' <4> user input of initial conditions\n\n');

    fprintf(' selection (1, 2, 3 or 4)\n\n');

    icflg = input('? ');

    if (icflg >= 1 && icflg <= 4)
        
        break;
        
    end
    
end

switch icflg
    
    case 1
        
        % periodic l1 orbit (tr 32-1168, pp 25,29; 74)

        y(1) = 0.300000161;
        y(3) = 0.0;
        y(2) = 0.0;
        y(4) = -2.536145497;

        ti = 0.0;
        tf = 5.349501906;

        mu = 0.012155092;

        % set plot boundaries

        xmin = -1.5;
        xmax = +1.5;
        ymin = -1.5;
        ymax = +1.5;

    case 2
        
        % periodic l2 orbit (tr 32-1168, pp 31,34; 126)

        y(1) = 2.840829343;
        y(3) = 0.0;
        y(2) = 0.0;
        y(4) = -2.747640074;

        ti = 0.0;
        tf = 2.0 * 5.966659294;

        mu = 0.012155085;

        % set plot boundaries

        xmin = -3.0;
        xmax = +3.0;
        ymin = -3.0;
        ymax = +3.0;

    case 3
        
        % periodic l3 orbit (tr 32-1168, pp 37,39; 63)

        y(1) = -1.600000312;
        y(3) = 0.0;
        y(2) = 0.0;
        y(4) = 2.066174572;

        ti = 0.0;
        tf = 2.0 * 3.151928156;

        mu = 0.012155092;

        % set plot boundaries

        xmin = -2.0;
        xmax = +2.0;
        ymin = -2.0;
        ymax = +2.0;

    case 4
        
        % user input of initial conditions

        fprintf('\nplease input the x component of the radius vector\n');

        y(1) = input('? ');

        fprintf('\nplease input the y component of the radius vector\n');

        y(3) = input('? ');

        fprintf('\nplease input the x component of the velocity vector\n');

        y(2) = input('? ');

        fprintf('\nplease input the y component of the velocity vector\n');

        y(4) = input('? ');

        while(1)

            fprintf('\nplease input the value for the earth-moon mass ratio\n');

            mu = input('? ');

            if (mu > 0.0)
                
                break;
                
            end
            
        end

        while(1)

            fprintf('\nplease input the normalized final time\n');

            tf = input('? ');

            if (abs(tf) > 0.0)
                
                break;
                
            end
        end

        % default values for plot boundaries

        xmin = -2.0;
        xmax = +2.0;
        ymin = -2.0;
        ymax = +2.0;

end

% request the integration step size

while(1)
    
    fprintf('\nplease input the normalized integration step size\n');

    fprintf('(a value of 0.01 is recommended)\n');

    dt = input('? ');

    if (abs(dt) > 0.0)
        
        break;
        
    end
    
end

% set ode45 options

options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);

% initialize

t2 = -dt;

npts = fix(tf / dt) + 1;

xplot = zeros(npts, 1);

yplot = zeros(npts, 1);

for i = 1:1:npts

    t1 = t2;

    t2 = t1 + dt;

    [twrk, ysol] = ode45(@crtbp_eqm, [t1, t2], y, options);

    xplot(i) = ysol(length(twrk), 1);

    yplot(i) = ysol(length(twrk), 3);

    y = ysol(length(twrk), 1:4);
    
end

% plot trajectory

plot(xplot, yplot, '-k', 'LineWidth', 1.5);

axis ([xmin xmax ymin ymax]);

axis square;

ylabel('normalized y coordinate');

xlabel('normalized x coordinate');

% label locations of Earth and Moon

hold on;

plot(-mu, 0, 'og', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

plot(1.0 - mu, 0, 'ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b');

% label libration points

switch icflg
    
    case 1
        
        plot(0.836892919, 0, '*r', 'MarkerSize', 8);

        title('Periodic Orbit about the L1 Libration Point', 'FontSize', 16);
        
    case 2
        
        plot(1.115699521, 0, '*r', 'MarkerSize', 8);

        title('Periodic Orbit about the L2 Libration Point', 'FontSize', 16);
        
    case 3
        
        plot(-1.005064527, 0, '*r', 'MarkerSize', 8);

        title('Periodic Orbit about the L3 Libration Point', 'FontSize', 16);
        
    case 4
        
        title('User Defined Initial Conditions', 'FontSize', 16);
end

% create tiff graphics disk file

print ('g3body.tif', '-dtiff'); 

fprintf('\n\n');


