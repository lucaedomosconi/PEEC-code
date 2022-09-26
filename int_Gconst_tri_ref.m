% Based on:
%
% Carley, Michael J.
% "Analytical formulae for potential integrals on triangles."
% Journal of Applied Mechanics 80, no. 4 (2013).

function I = int_Gconst_tri_ref(z, r1, r2, Theta)
  a = (r2 .* cos(Theta) - r1) ./ (r2 .* sin(Theta));
  phi = atan(a);
  beta2 = (r1.^2 + z.^2 .* (1 + a.^2)) ./ (1 + a.^2);
  alpha2 = z.^2 ./ beta2;
  alpha = sqrt(alpha2);
  beta = sqrt(beta2);
  I = beta .* (I0m11(Theta + phi, alpha) - I0m11(phi, alpha)) ...
      - abs(z) .* Theta;
end
