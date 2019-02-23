%
% USER-EDIT VARIABLES
%

% bearings run anticlockwise like in a unit circle

% runner start position
runnerPosition(1) = 100; % x
runnerPosition(2) = 100; % y
runnerPosition(3) = 100; % z

% runner velocity
runnerSpeed = 1;
% angle is converted to radians
runnerBearing = deg2rad(130); % XY
runnerElevation = deg2rad(0); % YZ

% chaser start position
chaserPosition(1) = 50; % x
chaserPosition(2) = 50; % y
chaserPosition(3) = 100; % z

% chaser speed
chaserSpeed = 2;

%
% CALCULATIONS
%

distanceVector = [chaserPosition(1) - runnerPosition(1) ...
    , chaserPosition(2) - runnerPosition(2) ...
    , chaserPosition(3) - runnerPosition(3)];
distance = norm(distanceVector); % magnitude (distance scalar value)

runnerVelocity = [runnerSpeed * cos(runnerBearing) * cos(runnerElevation) ...
    , runnerSpeed * sin(runnerBearing) * cos(runnerElevation) ...
    , runnerSpeed * sin(runnerElevation)];

% find time of collision using cosine rule
[timeUntilCollision1, timeUntilCollision2] = solveQuadratic(chaserSpeed^2 - runnerSpeed^2, 2*(dot(runnerVelocity, distanceVector)), -(distance^2));

% find points of collision

% collision point 1
collision1Position(1) = NaN; % x
collision1Position(2) = NaN; % y
collision1Position(3) = NaN; % z
if ~isnan(timeUntilCollision1) && timeUntilCollision1 > 0 % collision 1 exists
    collision1Position = runnerPosition + runnerVelocity*timeUntilCollision1;
end

% collision point 2
collision2Position(1) = NaN; % x
collision2Position(2) = NaN; % y
collision2Position(3) = NaN; % z
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

    % chaser bearing
    chaserBearing = atan(chaserVelocity(2)/chaserVelocity(1));
    fprintf('Chaser bearing for this collision: %f\n', rad2deg(chaserBearing));
    
    % chaser elevation
    d = sqrt(chaserVelocity(1)^2 * chaserVelocity(2)^2); % movement in xy plane
    if d == 0
        % account for one of them being 0
        d = max(chaserVelocity(1), chaserVelocity(2));
    end
    chaserElevation = atan(chaserVelocity(3) / d);
    fprintf('Chaser elevation for this collision: %f\n', rad2deg(chaserElevation));

    %
    % GRAPH
    %

    % plot runner quiver
    quiver3(runnerPosition(1), runnerPosition(2), runnerPosition(3) ...
        , runnerVelocity(1) * timeUntilClosestCollision ...
        , runnerVelocity(2) * timeUntilClosestCollision ...
        , runnerVelocity(3) * timeUntilClosestCollision ...
        , 'AutoScale', 'off');
    hold on
    
    % plot chaser quiver
    quiver3(chaserPosition(1), chaserPosition(2), chaserPosition(3) ...
        , chaserVelocity(1) * timeUntilClosestCollision ...
        , chaserVelocity(2) * timeUntilClosestCollision ...
        , chaserVelocity(3) * timeUntilClosestCollision ...
        , 'AutoScale', 'off');
    hold off
    
    % labels
    legend("runner", "chaser");
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
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
