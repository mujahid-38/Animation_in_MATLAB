% This program displays a mass attached to a horizontal spring and damper 
% fixed at the left end. The mass is released at the compressed spring 
% position due to which it oscillates and due to damping it stops its 
% oscillation after some time. It is a multiparametric analysis in which 
% you will observe that the greater mass results in the faster decay of 
% oscillation amplitude.

tic; clc; clear; close all;

% Setting up parameters for a 2D animation of a mass-spring-damper system

% Initial y-coordinate of the mass in the plane
mass_y_pos = 0;

% Global variable for the y-coordinate of the horizontal fixed end
global fixed_end_y_pos;
fixed_end_y_pos = -0.18;  % This horizontal fixed end becomes a line in the 2D animation

% Global variables for the x-coordinates of the horizontal fixed end
global fixed_end_x_pos1;
global fixed_end_x_pos2;
fixed_end_x_pos1 = 0;     % Starting x-coordinate of the horizontal fixed end
fixed_end_x_pos2 = 8;     % Ending x-coordinate of the horizontal fixed end

% Global variable for the width of the spring
global spring_width;
spring_width = 0.08;  % Width of the spring used in the animation


%% Simulation

%Simulation time
t_sim = 15;
t=0:0.3:t_sim;

% System Constant Parameters
m = 100; k = 100; b = [20 40 60]; F = 10;

X0 = [-.5,0];     % Initial conditions [position velocity]

% Solving the system with 3 different values of damping coefficient
[time1,y1]=ode45(@(t,y) ode_function(t,y,m,b(1),k,F), t, X0);
[time2,y2]=ode45(@(t,y) ode_function(t,y,m,b(2),k,F), t, X0);
[time3,y3]=ode45(@(t,y) ode_function(t,y,m,b(3),k,F), t, X0);

% Calculating the accelaration of system
a1 = gradient(y1(:,2));
a2 = gradient(y2(:,2));
a3 = gradient(y3(:,2));

% Find the maximum absolute acceleration among a1, a2, and a3
max_a = max(max([abs(a1), abs(a2), abs(a3)]));

% Find the maximum absolute velocity among y1(:,2), y2(:,2), and y3(:,2)
max_v = max(max([y1(:,2),y2(:,2),y3(:,2)]));

% Find the maximum absolute displacement among y1(:,1), y2(:,1), and y3(:,1)
max_x = max(max([y1(:,1),y2(:,1),y3(:,1)]));


%% Animation

% Set the position of the current figure window
set(gcf,'Position',[50 50 1080 480])

% Create a VideoWriter object for generating a video file
v = VideoWriter('single_mass_spring_damper.mp4','MPEG-4');

% Set the quality of the video (compression level)
v.Quality   = 50;

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
        y = y2;      % y2 is selected as current series for plotting
        a = a2;      % a2 is selected as current series for plotting
        col = 'b';
    elseif n==3
        y = y3;     % y3 is selected as current series for plotting
        a = a3;     % a3 is selected as current series for plotting
        col = 'g';
    end


    % Animation 
    for i=1:length(t)
        mass_x_pos = 1.5 + y(i,1);     % Shifting The Displacement values for animation 

        figure(1);clf
        nexttile        % Plot of Animated system
        hold on
        % Plotting the fixed end
        plot([fixed_end_x_pos1 fixed_end_x_pos2], [fixed_end_y_pos fixed_end_y_pos], 'LineWidth',4,'color','k');

        % Plotting the Mass
        plot(mass_x_pos,mass_y_pos,'s','MarkerSize',50,'MarkerFaceColor','r');

        % Plotting the damper
        damp_plot(0,mass_x_pos,-0.08);

        % Plotting the spring
        spring_plot(0,mass_x_pos,0.1);

        xlabel('x');ylabel('y');
        axis([0 8 -0.5 0.5]);
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

        axis([0 t_sim -1.1*max_x 1.1*max_x]);
        xlabel('Time [s]'),ylabel('Postion [m]')
        legend('b = 20', 'b = 40', 'b = 60')
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
        
        axis([0 t_sim -1.1*max_v 1.1*max_v]);
        xlabel('Time [s]'),ylabel('Velocity [ms^-^1]')
        legend('k = 30', 'k = 50', 'k = 70')
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
        
        axis([0 t_sim -1.1*max_a 1.1*max_a]);
        xlabel('Time [s]'),ylabel('Accelaration [ms^-^2]')
        legend('k = 30', 'k = 50', 'k = 70')
        hold off

        % Refresh the figure window to update its display
        refresh
        
        % Capture the current frame of the figure
        frame = getframe(gcf);

        % Check if the current time value (t1(i)) is equal to a specific simulation time (t_sim)
        if time1(i) == t_sim
            % If true, save the captured frame as an image file named 'dcl.jpg'
            imwrite(frame.cdata,'dcl.jpg');
        end

        % Write the captured frame to a video file 
        writeVideo(v,frame);
    end
end

% Close the VideoWriter object
close(v);

% Print the program execution time
toc;

% ode_function: Ordinary Differential Equation (ODE) function for a simple mass-spring-damper system
%
% Inputs:
%   - ~: The tilde (~) indicates that the function does not use the time variable explicitly
%   - y: The state vector [position; velocity]
%   - m: Mass of the system
%   - b: Damping coefficient
%   - k: Spring constant
%   - F: External force applied to the system
%
% Output:
%   - dy: Column vector containing the derivatives of the state variables

function dy = ode_function(~, y, m, b, k, F)
    % Equation for the first derivative of position (dy(1) / dt = y(2))
    dy(1) = y(2);

    % Equation for the first derivative of velocity
    % (dy(2) / dt = (1/m)*(F - k*y(1) - b*y(2)))
    dy(2) = (1/m) * (F - k*y(1) - b*y(2));

    % Return the derivatives as a column vector
    dy = dy';
end


% spring_plot: Function to plot a spring
%
% Inputs:
%   - x1: x-coordinate of the upper end of the spring
%   - x2: x-coordinate of the lower end of the spring
%   - y_pos: y-coordinate of the fixed end of the spring
%
% Globals:
%   - fixed_end_x_pos1: Global variable for the x-coordinate of the left fixed end
%
% Output:
%   - Plot of a horizontal spring between x1 and x2 with the fixed end at x_pos

function spring_plot(x1,x2,y_pos)
    % Access global variables
    global fixed_end_x_pos1
    global spring_width

    % Adjust x1 based on the left fixed end x-coordinate
    if x1 ~= fixed_end_x_pos1
        x1 = x1 + 0.22;
    end

    % Adjust x2 to extend the spring
    x2 = x2 - 0.55;

    % Create points for the spring
    x = linspace(x1,x2,12);
    y = spring_width*[0 0 -1 1 -1 1 -1 1 -1 1 0 0] + y_pos;
    
    % Plot the vertical spring with specified line width
    plot(x,y, 'LineWidth', 2)
end


% damp_plot: Function to plot a damper
%
% Inputs:
%   - x1: x-coordinate of the upper end of the damper
%   - x2: x-coordinate of the lower end of the damper
%   - y_pos: y-coordinate of the fixed end of the damper
%
% Globals:
%   - fixed_end_x_pos1: Global variable for the x-coordinate of the left fixed end
%
% Output:
%   - Plot of a horizontal damper between x1 and x2 with the fixed end at x_pos

function damp_plot(x1,x2,y_pos)

    % Access global variables
    global fixed_end_x_pos1

    % Adjust x1 based on the left fixed end x-coordinate
    if x1 ~= fixed_end_x_pos1
        x1 = x1 + 0.22;
    end

    % Adjust x2 of the damper
    x2 = x2 - 0.55;

    x = linspace(x1,x2,11);

    % heights of damper parts
    const1 = 0.04;
    const2 = 0.06;

    %plotting of damper parts
    plot([x(1),x(6)], [y_pos,y_pos], 'b', 'LineWidth', 2);
    plot([x(6),x(6)], [y_pos-const1,y_pos+const1], 'b', 'LineWidth', 2);
    plot([x(4), x(8), x(8), x(4)], [y_pos+const2,y_pos+const2,y_pos-const2,y_pos-const2], 'b', 'LineWidth', 2);
    plot([x(8),x(11)], [y_pos,y_pos], 'b', 'LineWidth', 2);

end
