% This program displays a simple pendulum on which a force of 2 N is
% applied. It is a multiparametric analysis in which you will observe that 
% the greater mass results in the smaller oscillation amplitude.


clc
clear
close all
tic;


%% Simulation

% Simulation time
t_sim = 0:0.05:5;

% System constant parameters
b = 0.1;
m = [0.5 1 2];
g = 9.81;
l = 1.5;
U = 2;
L = 1.5;

% Initial conditions [position velocity]
x0 = [0 0];

% Solving the system with 3 different values of mass
[t1,y1]=ode45(@(t,y) ode_function(t,y,m(1),g,l,b,U), t_sim, x0);
[t2,y2]=ode45(@(t,y) ode_function(t,y,m(2),g,l,b,U), t_sim, x0);
[t3,y3]=ode45(@(t,y) ode_function(t,y,m(3),g,l,b,U), t_sim, x0);


% Calculating the accelaration of system
a1 = gradient(y1(:,2));
a2 = gradient(y2(:,2));
a3 = gradient(y3(:,2));

% Find the maximum absolute acceleration among a1, a2, and a3
max_a = max(max([abs(a1), abs(a2), abs(a3)]));

% Find the maximum absolute velocity among y1(:,2), y2(:,2), and y3(:,2)
max_v = max(max([abs(y1(:,2)),abs(y2(:,2)),abs(y3(:,2))]));

% Find the maximum absolute displacement among y1(:,1), y2(:,1), and y3(:,1)
max_x = max(max([abs(y1(:,1)),abs(y2(:,1)),abs(y3(:,1))]));

if max_v == 0
    max_v = 1;
end
if max_a == 0
    max_a = 1;
end

%% Animation

% Set the position of the current figure window
set(gcf,'Position',[50 50 1080 560])

% Create a VideoWriter object for generating a video file
v = VideoWriter('simple_pendulum.mp4','MPEG-4');

% Set the quality of the video (compression level)
v.Quality   = 50;

% Open the VideoWriter object for writing
open(v);

for n = 1:1:3    % number of iterations based on the count of varying parameter
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
    for i = 1:length(t1)
    
        figure(1);clf
        subplot(2,2,1)        % Plot of Animated system
        ax = gca;
        ax.Position(4) = 1.3*ax.Position(4);
        hold on
            % Plotting the fixed end
            plot([-1 1], [0 0], 'LineWidth',5,'Color','k'); 
           
            % Plotting the arrow
            plot([-0.5 -0.25 -0.5],[-1.25 -1.5 -1.75], 'LineWidth',2,'Color','k')
            plot([-1 -0.25], [-1.5 -1.5], 'LineWidth',2,'Color','k')

            % Adding Text "F = force in newtons"
            text(-1.75,-1.5,sprintf('F = %.0f N',U),"FontWeight","bold")

            % Plotting the string
            plot([0 L*sin(y(i,1))],[0 -1*L*cos(y(i,1))],LineWidth=3)

            % Plotting the mass
            plot(L*sin(y(i,1)),-1*L*cos(y(i,1)),'o','MarkerFacecolor','b','MarkerSize',20)
            
            axis equal
            axis([-2 2 -2 2])
        hold off
        
        subplot(2,2,2)       % Plot of Angular Position
        tt = 1:i;
        hold on
            % Plotting all three series with default line width
            plot(t1,y1(:,1),'r-');
            plot(t1,y2(:,1),'b-');
            plot(t1,y3(:,1),'g-');
    
            % Plotting the current series with thicker line width
            plot(t1(tt),y(tt,1),col,'LineWidth',2)
    
            % Plotting the current point as a circle
            plot(t1(i),y(i,1),'o','MarkerFacecolor',col,'MarkerSize',5)
    
            legend('0.5 kg mass', '1 kg mass', '2 kg mass', Location='southeast')
            axis([0 t1(end) -1.5*max_x 1.5*max_x]);
            xlabel('Time [s]'),ylabel('Angular Postion [rad]')
        hold off
        
        subplot(2,2,3);       % Plot of Velocity
        hold on
            % Plotting all three series with default line width
            plot(t1,y1(:,2),'r-');
            plot(t1,y2(:,2),'b-');
            plot(t1,y3(:,2),'g-');
    
            % Plotting the current series with thicker line width
            plot(t1(tt),y(tt,2),col,'LineWidth',2)
    
            % Plotting the current point as a circle
            plot(t1(i),y(i,2),'o','MarkerFacecolor',col,'MarkerSize',5)
    
            axis([0 t1(end) -1.5*max_v 1.5*max_v]);
            xlabel('Time [s]'),ylabel('Angular Velocity [rad/s]')
        hold off
    
        subplot(2,2,4);       % Plot of Accelaration
        hold on   
            % Plotting all three series with default line width 
            plot(t1,a1,'r-');
            plot(t1,a2,'b-');
            plot(t1,a3,'g-');
    
            % Plotting the current series with thicker line width
            plot(t1(tt),a(tt),col,'LineWidth',2)
    
            % Plotting the current point as a circle
            plot(t1(i),a(i),'o','MarkerFacecolor',col,'MarkerSize',5)
    
            axis([0 t1(end) -1.5*max_a 1.5*max_a]);
            xlabel('Time [s]'),ylabel('Angular Accelaration [rad/s^2]')
        hold off
    
        % Writing data to the VideoWriter
        refresh

        % Capture the current frame of the figure
        frame = getframe(gcf);

        % Check if the current time value (t1(i)) is equal to a specific simulation time
        if t1(i) == t_sim(end)
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

% ode_function: Ordinary Differential Equation (ODE) function for a pendulum system
%
% Inputs:
%   - ~: The tilde (~) indicates that the function does not use the time variable explicitly
%   - x: State vector [theta; omega] representing the angle and angular velocity
%   - m: Mass of the pendulum bob
%   - g: Acceleration due to gravity
%   - l: Length of the pendulum
%   - b: Damping coefficient
%   - U: External force applied to the pendulum
%
% Output:
%   - xdot: Column vector containing the derivatives of the state variables

function xdot = ode_function(~, x, m, g, l, b, U)
    % Initialize the output vector
    xdot = zeros(2, 1);

    % Equation for the first derivative of angle (xdot(1) / dt = x(2))
    xdot(1) = x(2);

    % Equation for the first derivative of angular velocity
    % (xdot(2) / dt = (-b/(m*l^2))*x(2) - (g/l)*sin(x(1)) + (1/(m*l^2))*U)
    xdot(2) = (-b / (m * l^2)) * x(2) - (g / l) * sin(x(1)) + (1 / (m * l^2)) * U;
end
