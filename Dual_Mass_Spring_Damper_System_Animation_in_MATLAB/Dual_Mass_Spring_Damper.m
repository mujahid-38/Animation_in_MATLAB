% This program displays two masses attached to a horizontal springs and dampers 
% fixed at theleft end. A constant force is applied a one mass. It is a multiparametric analysis in which 
% you will observe that the greater applied force results in the greater displacement

tic; clc; clear; close all;

% Setting up parameters for a 2D animation of a mass-spring-damper system

% Initial y-coordinate of the mass in the plane
mass1_y_pos = 0;
mass2_y_pos = 0;

% Global variable for the y-coordinate of the horizontal fixed end
global fixed_end_y_pos;
fixed_end_y_pos = -0.18;  % This horizontal fixed end becomes a line in the 2D animation

% Global variables for the x-coordinates of the horizontal fixed end
global fixed_end_x_pos1;
global fixed_end_x_pos2;
fixed_end_x_pos1 = 0;     % Starting x-coordinate of the horizontal fixed end
fixed_end_x_pos2 = 3;     % Ending x-coordinate of the horizontal fixed end

% Global variable for the width of the spring
global spring_width;
spring_width = 0.08;  % Width of the spring used in the animation

%% Simulation


%Simulation time
t_sim=30;
t=0:0.4:t_sim;

% System Constant Parameters
m1=30;
m2=20;
k1=70;
k2=70;
b1=10;
b2=10;
F = [5 10 15];

% Initial conditions [position1 velocity1 position2 velocity2]
X0 = [0,0,0,0];   

%Solving the system with 3 different values of applied force
[t1,y1]=ode45(@(t,y) ode_function(t,y,m1,m2,b1,b2,k1,k2,F(1)), t, X0);
[t2,y2]=ode45(@(t,y) ode_function(t,y,m1,m2,b1,b2,k1,k2,F(2)), t, X0);
[t3,y3]=ode45(@(t,y) ode_function(t,y,m1,m2,b1,b2,k1,k2,F(3)), t, X0);

%Calculating the accelaration of the system
m1a1 = gradient(y1(:,2));
m2a1 = gradient(y1(:,4));
m1a2 = gradient(y2(:,2));
m2a2 = gradient(y2(:,4));
m1a3 = gradient(y3(:,2));
m2a3 = gradient(y3(:,4));

% Find the maximum absolute acceleration among m1a1, m2a1, m1a2, m2a2, m1a3 and m2a3
max_a = max(max([abs(m1a1), abs(m1a2), abs(m1a3), abs(m2a1), abs(m2a2), abs(m2a3)]));

% Find the maximum absolute velocity among y1(:,2), y2(:,2), y3(:,2), y1(:,4), y2(:,4), and y3(:,4)
max_v = max(max([y1(:,2),y1(:,4),y2(:,2),y2(:,4),y3(:,2),y3(:,4)]));

% Find the maximum absolute displacement among y1(:,1), y2(:,1), y3(:,1), y1(:,3), y2(:,3), and y3(:,3)
max_x = max(max([y1(:,1),y1(:,3),y2(:,1),y2(:,3),y3(:,1),y3(:,3)]));


%% Animation 

% Set the position of the current figure window
set(gcf,'Position',[50 50 1080 480])

% Create a VideoWriter object for generating a video file
v = VideoWriter('double_mass_spring_damper.mp4','MPEG-4');

% Set the quality of the video (compression level)
v.Quality   = 50;

% Open the VideoWriter object for writing
open(v);


for n = 1:1:3       % number of iterations based on the count of varying parameter
    % Selecting the current values to be highlighted with pointer in the
    % graphs
    if n == 1
        y = y1;     % y1 is selected as current series for plotting
        a1 = m1a1;  % m1a1 is selected as current series for plotting
        a2 = m2a1;  % m2a1 is selected as current series for plotting
        col1 = 'r';
        col2 = 'c';
    elseif n==2
        y = y2;     % y2 is selected as current series for plotting
        a1 = m1a2;   % m1a2 is selected as current series for plotting
        a2 = m2a2;   % m2a2 is selected as current series for plotting
        col1 = 'b';
        col2 = 'k';
    elseif n==3
        y = y3;     % y3 is selected as current series for plotting
        a1 = m1a3;  % m1a3 is selected as current series for plotting
        a2 = m2a3;  % m2a3 is selected as current series for plotting
        col1 = 'g';
        col2 = 'm';
    end

    for i=1:length(t)
        mass1_x_pos = 1.5 + y(i,1);      % Shifting The Displacement values for animation 
        mass2_x_pos = 0.5 + y(i,3);      % Shifting The Displacement values for animation 
        
        figure(1);clf
        tiledlayout(2,2,"TileSpacing","tight")
        nexttile    % Plot of Animated system
        hold on
        
        % Plotting the fixed end
        plot([fixed_end_x_pos1 fixed_end_x_pos2], [fixed_end_y_pos fixed_end_y_pos], 'LineWidth',4,'color','k');
        
        % Plotting the Masses
        plot(mass1_x_pos,mass1_y_pos,'s','MarkerSize',50,'MarkerFaceColor',col1);
        plot(mass2_x_pos,mass2_y_pos,'s','MarkerSize',50,'MarkerFaceColor',col2);
        
        % Plotting the damper
        damp_plot(fixed_end_x_pos1, mass2_x_pos, -0.08);
        damp_plot(mass2_x_pos, mass1_x_pos, -0.08);
        spring_plot(fixed_end_x_pos1, mass2_x_pos, 0.1);
        spring_plot(mass2_x_pos, mass1_x_pos, 0.1);
    
        xlabel('x');ylabel('y');
        axis([0 3 -0.5 0.5]);
        hold off
    
        nexttile;       % Plot of Position
        tt = 1:i;
        hold on
        % Plotting all of the series with default line width
        plot(t,y1(:,1),'r-');
        plot(t,y1(:,3),'c-');
        plot(t,y2(:,1),'b-');
        plot(t,y2(:,3),'k-');
        plot(t,y3(:,1),'g-');
        plot(t,y3(:,3),'m-');

        legend('v1 at 10 N', 'v2 at 10 N', 'v1 at 15 N', 'v2 at 15 N','v1 at 20 N', 'v2 at 20 N')

        % Plotting the current series with thicker line width
        plot(t(tt),y(tt,1),col1,'LineWidth',2)
        plot(t(tt),y(tt,3),col2,'LineWidth',2)

        % Plotting the current points as a circle
        plot(t(i),y(i,1),'o','MarkerFacecolor',col1,'MarkerSize',5)
        plot(t(i),y(i,3),'o','MarkerFacecolor',col2,'MarkerSize',5)

        legend('x1 at 10 N', 'x2 at 10 N', 'x1 at 15 N', 'x2 at 15 N','x1 at 20 N', 'x2 at 20 N', Location='southeast')
        axis([0 t_sim -1.1*max_x 1.25*max_x]);
        xlabel('Time [s]'),ylabel('Postion [m]')
        hold off
        
        nexttile;       % Plot of Velocity

        hold on

        % Plotting all of the series with default line width
        plot(t,y1(:,2),'r-');
        plot(t,y1(:,4),'c-');
        plot(t,y2(:,2),'b-');
        plot(t,y2(:,4),'k-');
        plot(t,y3(:,2),'g-');
        plot(t,y3(:,4),'m-');

        % Plotting the current series with thicker line width
        plot(t(tt),y(tt,2),col1,'LineWidth',2)
        plot(t(tt),y(tt,4),col2,'LineWidth',2)

        % Plotting the current point as a circle
        plot(t(i),y(i,2),'o','MarkerFacecolor',col1,'MarkerSize',5)
        plot(t(i),y(i,4),'o','MarkerFacecolor',col2,'MarkerSize',5)

        legend('v1 at 10 N', 'v2 at 10 N', 'v1 at 15 N', 'v2 at 15 N','v1 at 20 N', 'v2 at 20 N', Location='southeast')
        axis([0 t_sim -1.1*max_v 1.25*max_v]);
        xlabel('Time [s]'),ylabel('Velocity [ms^-^1]')

        nexttile;       % Plot of Accelaration
        hold on

        % Plotting all three series with default line width
        plot(t,m1a1,'r-');
        plot(t,m2a1,'c-');
        plot(t,m1a2,'b-');
        plot(t,m2a2,'k-');
        plot(t,m1a3,'g-');
        plot(t,m2a3,'m-');

        % Plotting the current series with thicker line width
        plot(t(tt),a1(tt),col1,'LineWidth',2)
        plot(t(tt),a2(tt),col2,'LineWidth',2)

        % Plotting the current point as a circle
        plot(t(i),a1(i),'o','MarkerFacecolor',col1,'MarkerSize',5)
        plot(t(i),a2(i),'o','MarkerFacecolor',col2,'MarkerSize',5)

        legend('a1 at 10 N', 'a2 at 10 N', 'a1 at 15 N', 'a2 at 15 N','a1 at 20 N', 'a2 at 20 N', Location='southeast')
        axis([0 t_sim -1.1*max_a 1.25*max_a]);
        xlabel('Time [s]'),ylabel('Accelaration [ms^-^2]')
        hold off
        
        % Refresh the figure window to update its display
        refresh

        % Capture the current frame of the figure
        frame = getframe(gcf);

        % Check if the current time value (t1(i)) is equal to a specific simulation time (t_sim)
        if t1(i) == t_sim
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


% ode_function: Ordinary Differential Equation (ODE) function for a double mass-spring-damper system
%
% Inputs:
%   - ~: The tilde (~) indicates that the function does not use the time variable explicitly
%   - y: The state vector [position; velocity]
%   - m1: Mass 1 of the system
%   - m2: Mass 2 of the system
%   - b1: Damping coefficient 1
%   - b2: Damping coefficient 2
%   - k1: Spring constant 1
%   - k2: Spring constant 2
%   - F: External force applied to the system
%
% Output:
%   - dy: Column vector containing the derivatives of the state variables

function dy=ode_function(~,y,m1,m2,b1,b2,k1,k2,F)
    % Equation for the first derivative of position of mass 1 (dy(1) / dt = y(2))
    dy(1)=y(2);

    % Equation for the first derivative of position of mass 2 (dy(1) / dt = y(2))
    dy(3)=y(4);

    % Equation for the first derivative of velocity of mass 1
    % (dy(2) / dt = (1/m1)*(F-k1*y(1) - b1*y(2) + k1*y(3) + b1*y(4))
    dy(2)=(1/m1)*(F-k1*y(1) - b1*y(2) + k1*y(3) + b1*y(4));

    % Equation for the first derivative of velocity of mass 2
    % (dy(2) / dt = (1/m2)*(k1*y(1) + b1*y(2) - (k1+k2)*y(3) - (b1+b2)*y(4))
    dy(4)=(1/m2)*(k1*y(1) + b1*y(2) - (k1+k2)*y(3) - (b1+b2)*y(4));

    % Return the derivatives as a column vector
    dy=dy';
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
    x2 = x2 - 0.22;

    % Create points for the spring
    x = linspace(x1,x2,12);
    y = 0.08*[0 0 -1 1 -1 1 -1 1 -1 1 0 0] + y_pos;

    % Plot the vertical spring with specified line width
    plot(x,y, 'LineWidth', 2,'Color',[0 0.5 0.8])
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
    x2 = x2 - 0.22;

    x = linspace(x1,x2,11);

    % heights of damper parts
    const1 = 0.04;
    const2 = 0.06;

    %plotting of damper parts
    plot([x(1),x(6)], [y_pos,y_pos], 'LineWidth', 2,'Color',[0.5 0 0.8]);
    plot([x(6),x(6)], [y_pos-const1,y_pos+const1], 'LineWidth', 2,'Color',[0.5 0 0.8]);
    plot([x(4), x(8), x(8), x(4)], [y_pos+const2,y_pos+const2,y_pos-const2,y_pos-const2], 'LineWidth', 2,'Color',[0.5 0 0.8]);
    plot([x(8),x(11)], [y_pos,y_pos], 'LineWidth', 2,'Color',[0.5 0 0.8]);

end
