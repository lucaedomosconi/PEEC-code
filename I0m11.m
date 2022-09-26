% Based on:
%
% Carley, Michael J.
% "Analytical formulae for potential integrals on triangles."
% Journal of Applied Mechanics 80, no. 4 (2013).

function v = I0m11(theta, alpha)
  alphap = sqrt(1 - alpha.^2);
  delta = sqrt(1 - (alpha .* sin(theta)).^2);
  v = alphap / 2 ...
      .* log((delta + alphap .* sin(theta)) ./ (delta - alphap .* sin(theta))) ...
      + alpha .* asin(alpha .* sin(theta));
end