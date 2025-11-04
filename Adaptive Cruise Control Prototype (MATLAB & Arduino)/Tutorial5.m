clc
clear

% % % % % % % % Define the function
% % % % % % % f = @(x) x .* tan(2 .* x) - 10;
% % % % % % % 
% % % % % % % % Plot the function in the interval [0, 5]
% % % % % % % fplot(f, [0 5])
% % % % % % % grid on
% % % % % % % xlabel('x')
% % % % % % % ylabel('f(x)')
% % % % % % % title('Plot of the function f(x) = x tan(2x) - 10')
% % % % % % % 
% % % % % % % % Use fzero to find the roots, providing initial guesses based on the plot
% % % % % % % root1 = fzero(f, 0.5);
% % % % % % % root2 = fzero(f, 1.5);
% % % % % % % root3 = fzero(f, 2.5);
% % % % % % % root4 = fzero(f, 3.5);
% % % % % % % root5 = fzero(f, 4.5);
% % % % % % % 
% % % % % % % % Display the roots
% % % % % % % roots = [root1, root2, root3, root4, root5];
% % % % % % % disp('Roots of the equation x tan(2x) = 10 are:')
% % % % % % % disp(roots)
% % % % % % 
% % % % % % % Define the matrix A
% % % % % % A = [0  -7  6;
% % % % % %      5  -4  3;
% % % % % %      10 -1  9;
% % % % % %      15  1  0;
% % % % % %      20  2 -1];
% % % % % % 
% % % % % % % Extract columns
% % % % % % time = A(:, 1)
% % % % % % force_col2 = A(:, 2)
% % % % % % force_col3 = A(:, 3)
% % % % % % 
% % % % % % % Plot columns 2 and 3 versus column 1
% % % % % % figure;
% % % % % % plot(time, force_col2, '-o', 'DisplayName', 'Column 2');
% % % % % % hold;
% % % % % % plot(time, force_col3, '-x', 'DisplayName', 'Column 3');
% % % % % % 
% % % % % % 
% % % % % % % Add labels and title
% % % % % % xlabel('Time (seconds)');
% % % % % % ylabel('Force (newtons)');
% % % % % % title('Force vs Time');
% % % % % % legend('show');
% % % % % % grid on;
% % % % % % 
% % % % % 
% % % % % % A. Given values
% % % % % v = 20; % m/s
% % % % % A = 25; % degrees
% % % % % g = 9.81; % m/s^2
% % % % % 
% % % % % % Convert angle to radians
% % % % % A_rad = deg2rad(A);
% % % % % 
% % % % % % Time of flight
% % % % % t_flight = 2 * v * sin(A_rad) / g;
% % % % % 
% % % % % % Maximum height
% % % % % h_max = (v^2 * sin(A_rad)^2) / (2 * g);
% % % % % 
% % % % % % Range
% % % % % range = (v^2 * sin(2 * A_rad)) / g;
% % % % % 
% % % % % % Display results
% % % % % fprintf('Time of flight: %.2f seconds\n', t_flight);
% % % % % fprintf('Maximum height: %.2f meters\n', h_max);
% % % % % fprintf('Range: %.2f meters\n', range);
% % % % % 
% % % % % %B. % Time vector
% % % % % t = linspace(0, t_flight, 100);
% % % % % 
% % % % % % Height and distance
% % % % % h = v * t * sin(A_rad) - 0.5 * g * t.^2;
% % % % % x = v * t * cos(A_rad);
% % % % % 
% % % % % % Plot
% % % % % figure;
% % % % % plot(x, h, 'b-', 'LineWidth', 2);
% % % % % xlabel('Horizontal Distance (m)');
% % % % % ylabel('Height (m)');
% % % % % title('Trajectory of the ball');
% % % % % grid on;
% % % % % 
% % % % % %C.% Given values
% % % % % A = 45; % degrees
% % % % % A_rad = deg2rad(A);
% % % % % 
% % % % % % Velocities
% % % % % velocities = [20, 24, 28, 32, 36];
% % % % % 
% % % % % % Plot
% % % % % figure;
% % % % % hold on;
% % % % % for v = velocities
% % % % %     % Time of flight
% % % % %     t_flight = 2 * v * sin(A_rad) / g;
% % % % % 
% % % % %     % Time vector
% % % % %     t = linspace(0, t_flight, 100);
% % % % % 
% % % % %     % Height and distance
% % % % %     h = v * t * sin(A_rad) - 0.5 * g * t.^2;
% % % % %     x = v * t * cos(A_rad);
% % % % % 
% % % % %     % Plot
% % % % %     plot(x, h, 'LineWidth', 2, 'DisplayName', ['v = ', num2str(v), ' m/s']);
% % % % % end
% % % % % hold off;
% % % % % 
% % % % % xlabel('Horizontal Distance (m)');
% % % % % ylabel('Height (m)');
% % % % % title('Trajectories of the ball for A = 45°');
% % % % % legend('show');
% % % % % grid on;
% % % % 
% % % % % Define the data
% % % % time = [1, 2, 3, 4, 5, 6, 7, 8, 10];
% % % % speed = [1210, 1866, 2301, 2564, 2724, 2881, 2879, 2915, 3010];
% % % % 
% % % % % Define the model function
% % % % modelFun = @(params, t) params(1) * (1 - exp(params(2) * t));
% % % % 
% % % % % Define the error function
% % % % errorFun = @(params) sum((speed - modelFun(params, time)).^2);
% % % % 
% % % % % Initial guesses for the parameters [b, c]
% % % % initialGuess = [3000, -0.5];
% % % % 
% % % % % Use fminsearch to minimize the error function
% % % % params = fminsearch(errorFun, initialGuess);
% % % % 
% % % % % Extract the parameters
% % % % b = params(1);
% % % % c = params(2);
% % % % 
% % % % % Display the parameters
% % % % fprintf('The value of b is: %.4f\n', b);
% % % % fprintf('The value of c is: %.4f\n', c);
% % % % 
% % % % % Plot the data and the fitted curve
% % % % fittedSpeed = modelFun(params, time);
% % % % figure;
% % % % plot(time, speed, 'o', 'DisplayName', 'Data');
% % % % hold on;
% % % % plot(time, fittedSpeed, '-r', 'DisplayName', sprintf('Fit: s(t) = %.4f (1 - exp(%.4f t))', b, c));
% % % % xlabel('Time (sec)');
% % % % ylabel('Speed (rpm)');
% % % % legend('show');
% % % % title('Motor Speed vs. Time');
% % % % grid on;
% % % 
% % % % Define the range of theta
% % % theta = linspace(0, 2*pi, 1000);
% % % 
% % % % Calculate the corresponding r values
% % % r = 6 * cos(0.8 * theta).^2 + theta;
% % % 
% % % % Create the polar plot
% % % figure;
% % % polarplot(theta, r);
% % % title('Polar Plot of r = 6 cos^2(0.8\theta) + \theta');
% % % grid on;
% % 
% % % Define the equations of the ellipses
% % ellipse1 = @(x, y) x.^2 + y.^2/4 - 1;
% % ellipse2 = @(x, y) 0.5833*x.^2 - 0.2887*x.*y + 0.4167*y.^2 - 1;
% % 
% % % Create a figure
% % figure;
% % 
% % % Use fimplicit to plot the first ellipse
% % fimplicit(ellipse1, [-2 2 -2 2], 'r');
% % hold on;
% % 
% % % Use fimplicit to plot the second ellipse
% % fimplicit(ellipse2, [-2 2 -2 2], 'b');
% % hold off;
% % 
% % % Label the axes
% % xlabel('x');
% % ylabel('y');
% % title('Intersection of Two Ellipses');
% % 
% % % Use ginput to find the intersection points
% % disp('Please click on the four intersection points in the plot');
% % [x, y] = ginput(4);
% % 
% % % Display the intersection points
% % intersection_points = [x, y];
% % disp('Intersection points:');
% % disp(intersection_points);
% 
% % Define the parameters
% a = 1;  % Radius of the helix
% t = linspace(0, 10*pi, 1000);  % Parameter t from 0 to 10*pi
% 
% % Define the values of b for the three cases
% b_values = [0.1, 0.2, -0.1];
% 
% % Create a new figure
% figure;
% 
% % Loop over each value of b and plot the corresponding helix
% for i = 1:length(b_values)
%     b = b_values(i);
% 
%     % Parametric equations for the helix
%     x = a * cos(t);
%     y = a * sin(t);
%     z = b * t;
% 
%     % Create a subplot for each value of b
%     subplot(3, 1, i);
%     plot3(x, y, z, 'LineWidth', 2);
%     grid on;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     title(['Helix with b = ', num2str(b)]);
%     axis equal;
% end
% 
% % Adjust the layout
% sgtitle('Three-dimensional Plots of the Helical Paths');

% Define the temperature distribution function
T = @(x, y) 80 * exp(-(x - 1).^2) .* exp(-3 * (y - 1).^2);

% Create a grid for x and y
[x, y] = meshgrid(linspace(0, 2, 100), linspace(0, 2, 100));

% Calculate the temperature values over the grid
temperature = T(x, y);

% Plot the surface plot
figure;
surf(x, y, temperature);
xlabel('x');
ylabel('y');
zlabel('Temperature (°C)');
title('Surface Plot of Temperature Distribution');
colorbar;
shading interp;

% Plot the contour plot
figure;
contour(x, y, temperature, 20);
xlabel('x');
ylabel('y');
title('Contour Plot of Temperature Distribution');
colorbar;

% Evaluate the temperature at the corner x = y = 0
corner_temperature = T(0, 0);
disp(['The temperature at the corner (x, y) = (0, 0) is ', num2str(corner_temperature), ' °C']);
