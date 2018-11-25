%
% USER-EDIT VARIABLES
%

% runner start position
runnerPosition(1) = 100;
runnerPosition(2) = 100;

% runner velocity
runnerSpeed = 1;
runnerDirection = deg2rad(100); % angle is converted to radians

% chaser start position
chaserPosition(1) = 50;
chaserPosition(2) = 50;

% chaser speed
chaserSpeed = 2;

%
% CALCULATIONS
%

distanceVector = [chaserPosition(1) - runnerPosition(1), chaserPosition(2) - runnerPosition(2)];
distance = norm(distanceVector); % magnitude (distance scalar value)

runnerVelocity = [runnerSpeed * cos(runnerDirection), runnerSpeed * sin(runnerDirection)];

% find time of collision
[timeUntilCollision1, timeUntilCollision2] = solveQuadratic(chaserSpeed^2 - runnerSpeed^2, 2*(dot(runnerVelocity, distanceVector)), -(distance^2));

% find points of collision

% collision point 1
collision1Position(1) = NaN;
collision1Position(2) = NaN;
if ~isnan(timeUntilCollision1) && timeUntilCollision1 > 0
    collision1Position = runnerPosition + runnerVelocity*timeUntilCollision1;
end

% collision point 2
collision2Position(1) = NaN;
collision2Position(2) = NaN;
if ~isnan(timeUntilCollision2) && timeUntilCollision2 > 0
    collision2Position = runnerPosition + runnerVelocity*timeUntilCollision2;
end

% chaser velocity
chaserVelocity = NaN;
% base off of a valid collision
if ~isnan(collision1Position(1))
    chaserVelocity = (collision1Position - chaserPosition) / timeUntilCollision1;
elseif ~isnan(collision2Position(1))
    chaserVelocity = (collision2Position - chaserPosition) / timeUntilCollision2;
end

% chaser direction
if ~isnan(chaserVelocity(1))
    chaserDirection = atan(chaserVelocity(2)/chaserVelocity(1));
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
