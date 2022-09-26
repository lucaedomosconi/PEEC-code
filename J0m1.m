% Based on:
%
% Carley, Michael J.
% "Analytical formulae for potential integrals on triangles."
% Journal of Applied Mechanics 80, no. 4 (2013).

function v = J0m1(theta)
  v = log((1 + sin(theta)) ./ (1 - sin(theta)));
end
