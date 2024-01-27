function out = LK_RandomPointsInCircle(myR)
%
% LK_RandomPointsInCircle creates random points drawn from a circular area.
%
% Input: structure with fields
%   maxR    --> maximum radius
%   minR    --> minimum radius
%   N       --> number of points to create
%   centerX --> x-center of circle
%   centerY --> y-center of circle
%
% Output:
%   (x/y)-locations drawn from a circular area.
%
% Lukas Kunz, 2023

% limits in the range [a, 1]
limits  = ([myR.minR, myR.maxR] ./ myR.maxR) .^ 2; % cave: square

% calculate random points within circle in polar coordinates
r       = sqrt(rand(myR.N, 1) .* range(limits) + min(limits)); % cave: sqrt
theta   = rand(myR.N, 1) .* 2 .* pi;

% transform into desired range
rTransf = r .* myR.maxR;

% transform into Cartesian coordinates
x       = myR.centerX + rTransf .* cos(theta);
y       = myR.centerY + rTransf .* sin(theta);

% output
out     = [x, y];

end