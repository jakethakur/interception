%
% USER-EDIT VARIABLES
%

% bearings run anticlockwise like in a unit circle

% runner start position
runnerPosition(1) = 100; % x
runnerPosition(2) = 100; % y

% runner velocity
runnerSpeed = 1;
runnerDirection = deg2rad(130); % angle is converted to radians

% chaser start position
chaserPosition(1) = 50; % x
chaserPosition(2) = 50; % y

% chaser speed
chaserSpeed = 2;

%
% CALCULATIONS
%

distanceVector = [chaserPosition(1) - runnerPosition(1), chaserPosition(2) - runnerPosition(2)];
distance = norm(distanceVector); % magnitude (distance scalar value)

runnerVelocity = [runnerSpeed * cos(runnerDirection), runnerSpeed * sin(runnerDirection)];

% find time of collision using cosine rule
[timeUntilCollision1, timeUntilCollision2] = solveQuadratic(chaserSpeed^2 - runnerSpeed^2, 2*(dot(runnerVelocity, distanceVector)), -(distance^2));

% find points of collision

% collision point 1
collision1Position(1) = NaN; % x
collision1Position(2) = NaN; % y
if ~isnan(timeUntilCollision1) && timeUntilCollision1 > 0 % collision 1 exists
    collision1Position = runnerPosition + runnerVelocity*timeUntilCollision1;
end

% collision point 2
collision2Position(1) = NaN; % x
collision2Position(2) = NaN; % y
if ~isnan(timeUntilCollision2) && timeUntilCollision2 > 0 % collision 2 exists
    collision2Position = runnerPosition + runnerVelocity*timeUntilCollision2;
end

% find the closest valid collision
closestCollisionPosition = NaN;
timeUntilClosestCollision = NaN;
if isnan(collision1Position(1)) && isnan(collision2Position(1)) % neither exist
    
elseif isnan(collision2Position(1)) && ~isnan(collision1Position(1)) % just collision 1 exists
    closestCollisionPosition = collision1Position;
    timeUntilClosestCollision = timeUntilCollision1;
    
elseif isnan(collision1Position(1)) && ~isnan(collision2Position(1)) % just collision 2 exists
    closestCollisionPosition = collision2Position;
    timeUntilClosestCollision = timeUntilCollision2;
    
elseif timeUntilCollision1 < timeUntilCollision2 % collision 1 is before collision 2
    closestCollisionPosition = collision1Position;
    timeUntilClosestCollision = timeUntilCollision1;
    
elseif timeUntilCollision2 < timeUntilCollision1 % collision 2 is before collision 1
    closestCollisionPosition = collision2Position;
    timeUntilClosestCollision = timeUntilCollision2;
    
end

if ~isnan(closestCollisionPosition)
    % there has been a valid collision
    
    % print coords
    % https://stackoverflow.com/a/27841544/9713957
    g = sprintf('%f ', closestCollisionPosition); % convert vector to string first
    fprintf('A valid collision has been found: %s\n', g);

    % chaser velocity (base off of the closest valid collision)
    chaserVelocity = (closestCollisionPosition - chaserPosition) / timeUntilClosestCollision;

    % chaser direction
    chaserDirection = atan(chaserVelocity(2)/chaserVelocity(1));
    fprintf('Chaser direction for this collision: %f\n', rad2deg(chaserDirection));

    %
    % GRAPH
    %

    % plot runner quiver
    quiver(runnerPosition(1), runnerPosition(2), runnerVelocity(1) * timeUntilClosestCollision, runnerVelocity(2) * timeUntilClosestCollision, 'AutoScale', 'off');
    hold on
    
    % plot chaser quiver
    quiver(chaserPosition(1), chaserPosition(2), chaserVelocity(1) * timeUntilClosestCollision, chaserVelocity(2) * timeUntilClosestCollision, 'AutoScale', 'off');
    hold off
    
    legend("runner", "chaser");
    
else
    % no valid collision has been found
    disp("No valid collision has been found.");
end

%
% ANIMATION
%

% https://uk.mathworks.com/help/matlab/creating_plots/trace-marker-along-line.html
% https://uk.mathworks.com/help/matlab/creating_plots/move-group-of-objects-along-line.html

%
% FUNCTIONS
%
function [root1, root2] = solveQuadratic(a, b, c)

  d = b^2 - 4*a*c; % your number under the root sign in quad. formula

  % real numbered distinct roots?
  if d > 0
    root1 = (-b+sqrt(d))/(2*a);
    root2 = (-b-sqrt(d))/(2*a);
  % real numbered degenerate root?
  elseif d == 0 
    root1 = -b/(2*a);
    root2 = NaN;
  % complex roots, return NaN, NaN
  else
    root1 = NaN;
    root2 = NaN;
  end    
end
