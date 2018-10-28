%
% SIMULATION VARIABLES
%

simulationSpeed = 1; % speed of objects is multiplied by this (CURRENTLY DOES NOTHING)
simulationLength = 1000; % number of iterations simulation should be run for

%
% SETUP
%

runnerPositionX = zeros(simulationLength, 1);
runnerPositionY = zeros(simulationLength, 1);

chaserPositionX = zeros(simulationLength, 1);
chaserPositionY = zeros(simulationLength, 1);

collision = false; % if there has been a collision
collisionX;
collisionY;

%
% USER-EDIT VARIABLES
%

% runner start position
runnerPositionX(1, 1) = 100;
runnerPositionY(1, 1) = 100;

% runner velocity
runnerSpeed = 1;
runnerDirection = deg2rad(100); % angle is converted to radians

% chaser start position
chaserPositionX(1, 1) = 200;
chaserPositionY(1, 1) = 200;

% chaser speed
chaserSpeed = 2;

%
% DO NOT EDIT THESE VARIABLES
%

% chaser direction
chaserDirection = deg2rad(200); % angle is converted to radians

%
% CALCULATIONS
%

% iterate for simulation length
for i = 2:simulationLength+1
   % move runner
   runnerPositionX(i, 1) = runnerPositionX(i-1, 1) + cos(runnerDirection) * runnerSpeed;
   runnerPositionY(i, 1) = runnerPositionY(i-1, 1) + sin(runnerDirection) * runnerSpeed;
   
   % move chaser
   chaserPositionX(i, 1) = chaserPositionX(i-1, 1) + cos(chaserDirection) * chaserSpeed;
   chaserPositionY(i, 1) = chaserPositionY(i-1, 1) + sin(chaserDirection) * chaserSpeed;
   
   % check for collision
   if chaserPositionX(i, 1) == runnerPositionX(i, 1) && chaserPositionY(i, 1) == runnerPositionY(i, 1)
       collision = true;
       collisionX = chaserPositionX(i, 1);
       collisionY = chaserPositionY(i, 1);
   end
end

%
% PLOT
%

figure;
plot(runnerPositionX, runnerPositionY, chaserPositionX, chaserPositionY);

%
% ANIMATION
%

% https://uk.mathworks.com/help/matlab/creating_plots/trace-marker-along-line.html
% https://uk.mathworks.com/help/matlab/creating_plots/move-group-of-objects-along-line.html