function [p, w] = tri_gauss_points(n, X1, X2, X3)
%   This function evaluates \iint_K f(x,y) dxdy using
%   the Gaussian quadrature of order N where K is a
%   triangle with vertices X1, X2 and X3.
  [eta, v] = tri_gauss_points_ref(n);  % get quadrature points and weights
  eta1 = reshape(eta(:, 1), [1, 1, size(eta, 1)]);
  eta2 = reshape(eta(:, 2), [1, 1, size(eta, 1)]);
  v = reshape(v, [1, 1, size(v, 1)]);
  % calculate the area of the triangle
  A = norm(cross(X2-X1, X3-X1, 1), 2, 'cols') / 2;
  % Gauss points and weights
  p = X1.*(1 - eta1 - eta2) + X2.*eta1 + X3.*eta2;
  w = v .* A;
end