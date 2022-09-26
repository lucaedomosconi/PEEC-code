% Based on:
%
% Carley, Michael J.
% "Analytical formulae for potential integrals on triangles."
% Journal of Applied Mechanics 80, no. 4 (2013).

function I = int_Gconst_tri_mesh(nodes, connectivity, p)
  n = size(connectivity, 2);
  elem = nodes(:,connectivity(:));
  elem = reshape(elem, [3, 3, n]);
  v1 = squeeze(elem(:,2,:) - elem(:,1,:));
  v2 = squeeze(elem(:,3,:) - elem(:,1,:));
  n = cross(v1, v2, 1);
  n ./= norm(n, 2, 'cols');
  clear v1
  clear v2
  z = dot(squeeze(p - elem(:,1,:)), n);
  origin = p - n .* repmat(z, 3, 1);
  l1 = origin - squeeze(elem(:,1,:));
  l2 = origin - squeeze(elem(:,2,:));
  l3 = origin - squeeze(elem(:,3,:));
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
  I(idx1)  = int_Gconst_tri_ref(z(idx1), r1(idx1), r2(idx1), theta1(idx1));
  I(idx2) += int_Gconst_tri_ref(z(idx2), r2(idx2), r3(idx2), theta2(idx2));
  I(idx3) += int_Gconst_tri_ref(z(idx3), r3(idx3), r1(idx3), theta3(idx3));
end