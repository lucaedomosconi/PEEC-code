% Based on:
%
% Carley, Michael J.
% "Analytical formulae for potential integrals on triangles."
% Journal of Applied Mechanics 80, no. 4 (2013).

function I = int_Gconst_tri(X1, X2, X3, p)
  n = cross(X2 - X1, X3 - X1, 1);
  n ./= norm(n, 2);
  n = repmat(n, 1, size(p, 2));
  z = dot(p - X1, n);
  origin = p - z .* n;
  l1 = origin - X1;
  l2 = origin - X2;
  l3 = origin - X3;
  r1 = norm(l1, 2, 'cols');
  r2 = norm(l2, 2, 'cols');
  r3 = norm(l3, 2, 'cols');
  theta1 = dot(cross(l1, l2, 1), n, 1);
  theta2 = dot(cross(l2, l3, 1), n, 1);
  theta3 = dot(cross(l3, l1, 1), n, 1);
  idx1 = find(abs(theta1) > 1e-10 & ~isnan(theta1));
  idx2 = find(abs(theta2) > 1e-10 & ~isnan(theta2));
  idx3 = find(abs(theta3) > 1e-10 & ~isnan(theta3));
  theta1 = sign(theta1);
  theta2 = sign(theta2);
  theta3 = sign(theta3);  
  theta1 .*= acos(dot(l1, l2) ./ (r1 .* r2));
  theta2 .*= acos(dot(l2, l3) ./ (r2 .* r3));
  theta3 .*= acos(dot(l3, l1) ./ (r3 .* r1));
  I = zeros(size(theta1));
  I(idx1)  = int_Gconst_tri_ref(z(idx1), r1(idx1), r2(idx1), theta1(idx1));
  I(idx2) += int_Gconst_tri_ref(z(idx2), r2(idx2), r3(idx2), theta2(idx2));
  I(idx3) += int_Gconst_tri_ref(z(idx3), r3(idx3), r1(idx3), theta3(idx3));
end