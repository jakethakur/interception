%
% INPUT VARIABLES
%

% bearings run anticlockwise like in a unit circle

% target start position
targetPosition(1) = 100; % x
targetPosition(2) = 100; % y
targetPosition(3) = 0; % z

% target velocity
targetSpeed = 1;
% angle is converted to radians
targetBearing = deg2rad(130); % XY
targetElevation = deg2rad(100); % YZ

% interceptor start position
interceptorPosition(1) = 50; % x
interceptorPosition(2) = 50; % y
interceptorPosition(3) = 100; % z

% interceptor speed
interceptorSpeed = 2;

%
% CALCULATIONS
%

distanceVector = [interceptorPosition(1) - targetPosition(1) ...
    , interceptorPosition(2) - targetPosition(2) ...
    , interceptorPosition(3) - targetPosition(3)];
distance = norm(distanceVector); % magnitude (distance scalar value)

targetVelocity = [targetSpeed * cos(targetBearing) * cos(targetElevation) ...
    , targetSpeed * sin(targetBearing) * cos(targetElevation) ...
    , targetSpeed * sin(targetElevation)];

% find time of collision using cosine rule
[timeUntilCollision1, timeUntilCollision2] = solveQuadratic(interceptorSpeed^2 - targetSpeed^2, 2*(dot(targetVelocity, distanceVector)), -(distance^2));

% find points of collision

% collision point 1
collision1Position(1) = NaN; % x
collision1Position(2) = NaN; % y
collision1Position(3) = NaN; % z
if ~isnan(timeUntilCollision1) && timeUntilCollision1 > 0 % collision 1 exists
    collision1Position = targetPosition + targetVelocity*timeUntilCollision1;
end

% collision point 2
collision2Position(1) = NaN; % x
collision2Position(2) = NaN; % y
collision2Position(3) = NaN; % z
if ~isnan(timeUntilCollision2) && timeUntilCollision2 > 0 % collision 2 exists
    collision2Position = targetPosition + targetVelocity*timeUntilCollision2;
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
    
    % print time taken for collisoin
    fprintf('Time taken for collision: %f\n', timeUntilClosestCollision);

    % interceptor velocity (base off of the closest valid collision)
    interceptorVelocity = (closestCollisionPosition - interceptorPosition) / timeUntilClosestCollision;

    % interceptor bearing
    interceptorBearing = atan(interceptorVelocity(2)/interceptorVelocity(1));
    fprintf('Interceptor bearing for this collision: %f\n', rad2deg(interceptorBearing));
    
    % interceptor elevation
    d = sqrt(interceptorVelocity(1)^2 * interceptorVelocity(2)^2); % movement in xy plane
    if d == 0
        % account for one of them being 0
        d = max(interceptorVelocity(1), interceptorVelocity(2));
    end
    interceptorElevation = atan(interceptorVelocity(3) / d);
    fprintf('Interceptor elevation for this collision: %f\n', rad2deg(interceptorElevation));

    %
    % GRAPH
    %

    % plot target quiver
    quiver3(targetPosition(1), targetPosition(2), targetPosition(3) ...
        , targetVelocity(1) * timeUntilClosestCollision ...
        , targetVelocity(2) * timeUntilClosestCollision ...
        , targetVelocity(3) * timeUntilClosestCollision ...
        , 'AutoScale', 'off');
    hold on
    
    % plot interceptor quiver
    quiver3(interceptorPosition(1), interceptorPosition(2), interceptorPosition(3) ...
        , interceptorVelocity(1) * timeUntilClosestCollision ...
        , interceptorVelocity(2) * timeUntilClosestCollision ...
        , interceptorVelocity(3) * timeUntilClosestCollision ...
        , 'AutoScale', 'off');
    hold off
    
    % labels
    legend("target", "interceptor");
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
else
    % no valid collision has been found
    disp("No valid collision has been found.");
end

%
% FUNCTIONS
%

% https://stackoverflow.com/a/34828707/9713957
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
