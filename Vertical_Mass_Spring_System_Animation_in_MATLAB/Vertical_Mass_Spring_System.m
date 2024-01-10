% This program displays a mass attached to a vetical spring fixed at the
% upper end. The mass is released at the position higher than the
% unstretched position (compressed state) and it goes downwards due to 
% the force of gravity. It is a multiparametric analysis in which you will
% observe that by increasing the hung mass, its steady state position is
% lowered because of its weight.

tic; clc; clear; close all;

%Simulation time
t_sim = 15;
           
mass_x_pos = 0;     % x-coordinate of mass in the plane
global fixed_end_y_pos;
fixed_end_y_pos = 10;       % y-coordinate of horizontal fixed end (becomes a line in 2D animation)
fixed_end_x_pos1 = -5;      % x-coordinate of starting point of horizontal fixed end
fixed_end_x_pos2 = 5;       % x-coordinate of ending point of horizontal fixed end

global spring_width;
spring_width = 0.4;

% System Constant Parameters
m = [10 20 30]; k = 200; 
b = 10; F = 50;

%% Simulation             
t=0:0.05:t_sim;
              
              
X0 = [3,0];     % Initial conditions [position velocity]

%Solving the system with 3 different values of mass
[t1,y1]=ode45(@(t,y) ode_function(t,y,m(1),b,k,F), t, X0);
[t2,y2]=ode45(@(t,y) ode_function(t,y,m(2),b,k,F), t, X0);
[t3,y3]=ode45(@(t,y) ode_function(t,y,m(3),b,k,F), t, X0);

% Calculating the accelaration of system
a1 = gradient(y1(:,2));
a2 = gradient(y2(:,2));
a3 = gradient(y3(:,2));


% Find the maximum absolute acceleration among a1, a2, and a3
max_a = max(max([abs(a1), abs(a2), abs(a3)]));

% Find the maximum absolute velocity among y1(:,2), y2(:,2), and y3(:,2)
max_v = max(max([abs(y1(:,2)), abs(y2(:,2)), abs(y3(:,2))]));

% Find the maximum absolute displacement among y1(:,1), y2(:,1), and y3(:,1)
max_x = max(max([abs(y1(:,1)), abs(y2(:,1)), abs(y3(:,1))]));

%% Animation

% Set the position of the current figure window
set(gcf, 'Position', [50 50 1080 480])

% Create a VideoWriter object for generating a video file
v = VideoWriter('vertical_mass_spring_system_1.mp4', 'MPEG-4');

% Set the quality of the video (compression level)
v.Quality = 50;

% Open the VideoWriter object for writing
open(v);

for n = 1:1:3   % number of iterations based on the count of varying parameter
    % Selecting the current values to be highlighted with pointer in the
    % graphs
    if n == 1   
        y = y1;     % y1 is selected as current series for plotting
        a = a1;     % a1 is selected as current series for plotting
        col = 'r';
    elseif n==2
        y = y2;     % y2 is selected as current series for plotting
        a = a2;     % a2 is selected as current series for plotting
        col = 'b';
    elseif n==3
        y = y3;     % y3 is selected as current series for plotting
        a = a3;     % a3 is selected as current series for plotting
        col = 'g';
    end


    % Animation 
    for i=1:length(t)
        mass_y_pos = 1.5 + y(i,1);     % Shifting The Displacement values for animation 
    
        figure(1);clf
        nexttile        % Plot of Animated system
        hold on
        % Plotting the fixed end
        plot([fixed_end_x_pos1 fixed_end_x_pos2], [fixed_end_y_pos fixed_end_y_pos], 'LineWidth', 4, 'color', 'k');
        
        % Plotting the Mass
        plot(mass_x_pos, mass_y_pos, 's', 'MarkerSize', 40, 'MarkerFaceColor', 'b');
        
        % Plotting the spring
        vertical_spring_plot(fixed_end_y_pos,mass_y_pos,mass_x_pos)

        xlabel('x');ylabel('y');
        axis([-5 5 -10 10]);
        hold off
    
        nexttile;       % Plot of Position
        tt = 1:i;
        hold on
        % Plotting all three series with default line width
        plot(t,y1(:,1),'r-');
        plot(t,y2(:,1),'b-');
        plot(t,y3(:,1),'g-');

        % Plotting the current series with thicker line width
        plot(t(tt),y(tt,1),col,'LineWidth',2)

        % Plotting the current point as a circle
        plot(t(i),y(i,1),'o','MarkerFacecolor',col,'MarkerSize',5)

        axis([0 t_sim -max_x-1 max_x+1]);
        xlabel('Time [s]'),ylabel('Postion [m]')
        legend('m = 10 kg', 'm = 20 kg', 'm = 30 kg')
        hold off
        
        nexttile;       % Plot of Velocity
        hold on
        
        % Plotting all three series with default line width
        plot(t,y1(:,2),'r-');
        plot(t,y2(:,2),'b-');
        plot(t,y3(:,2),'g-');

        % Plotting the current series with thicker line width
        plot(t(tt),y(tt,2),col,'LineWidth',2)

        % Plotting the current point as a circle
        plot(t(i),y(i,2),'o','MarkerFacecolor',col,'MarkerSize',5)

        axis([0 t_sim -max_v-1 max_v+1]);
        xlabel('Time [s]'),ylabel('Velocity [ms^-^1]')
        hold off
    
        nexttile;       % Plot of Accelaration
        hold on  

        % Plotting all three series with default line width  
        plot(t,a1,'r-');
        plot(t,a2,'b-');
        plot(t,a3,'g-');

        % Plotting the current series with thicker line width
        plot(t(tt),a(tt),col,'LineWidth',2)

        % Plotting the current point as a circle
        plot(t(i),a(i),'o','MarkerFacecolor',col,'MarkerSize',5)

        axis([0 t_sim -max_a-1 max_a+1]);
        xlabel('Time [s]'),ylabel('Accelaration [ms^-^2]')
        hold off

        % Refresh the figure window to update its display
        refresh;

        % Capture the current frame of the figure
        frame = getframe(gcf);

        % Check if the current time value (t1(i)) is equal to a specific simulation time (t_sim)
        if t1(i) == t_sim
            % If true, save the captured frame as an image file named 'dcl.jpg'
            imwrite(frame.cdata, 'dcl.jpg');
        end

        % Write the captured frame to a video file 
        writeVideo(v, frame);
    end
end
% Close the VideoWriter object
close(v);

% Print the time duration of the program execution
toc;


% Define an ODE function for a mass-spring-damper system
% Inputs:
%   - ~: The tilde (~) indicates that the function does not use the time variable explicitly
%   - y: The state vector [position; velocity]
%   - m: Mass of the system
%   - b: Damping coefficient
%   - k: Spring constant
%   - F: External force applied to the system

function dy = ode_function(~, y, m, b, k, F)
    % Equation for the first derivative of position (dy(1) / dt = y(2))
    dy(1) = y(2);

    % Equation for the first derivative of velocity
    % (dy(2) / dt = (1/m)*(-F - m*g - k*y(1) - b*y(2)))
    dy(2) = (1/m) * (-F - m*9.81 - k*y(1) - b*y(2));

    % Return the derivatives as a column vector
    dy = dy';
end


% vertical_spring_plot: Function to plot a vertical spring
%
% Inputs:
%   - y1: y-coordinate of the upper end of the spring
%   - y2: y-coordinate of the lower end of the spring
%   - x_pos: x-coordinate of the fixed end of the spring
%
% Globals:
%   - fixed_end_y_pos: Global variable for the y-coordinate of the fixed end
%   - spring_width: Global variable for the width of the spring
%
% Output:
%   - Plot of a vertical spring between y1 and y2 with the fixed end at x_pos

function vertical_spring_plot(y1, y2, x_pos)
    % Access global variables
    global fixed_end_y_pos
    global spring_width

    % Adjust y1 based on the fixed end y-coordinate
    if y1 ~= fixed_end_y_pos
        y1 = y1 - 2.6;
    else
        y1 = y1 - 0.2;
    end

    % Adjust y2 to extend the spring
    y2 = y2 + 2.6;

    % Create points for the spring
    y = linspace(y1, y2, 12);
    x = spring_width * [0 0 -1 1 -1 1 -1 1 -1 1 0 0] + x_pos;

    % Plot the vertical spring with specified line width
    plot(x, y, 'LineWidth', 2);
end
