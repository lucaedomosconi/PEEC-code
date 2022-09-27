clearvars
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 4.6729e6;
freq = 1000;
thickness = 3e-3;
omega = 2 * pi * freq;
mu0 = 4 * pi * 1e-7;
mr = 1.;
mu = mu0*mr;
skindepth = sqrt(2/omega/mu/sigma);

% uniform background magnetic field (peak values)
B0x = 0;
B0y = 0;
B0z = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh construction                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0.1;
b = 0.1;
m = 20;
n = 20;
figure
[Y,X] = meshgrid(linspace(-b/2,b/2,n+1), linspace(-a/2,a/2,m+1));
%Z = zeros(size(X));
Z = 0. * (X.^2 - Y.^2);
% nodes matrix: (size 3 x num_nodes)
nodes = [X(:)'; Y(:)'; Z(:)'];

% connectivity matrix: (size 3 x num_elements)
cnc = [1,   m+3;
       2,   m+2;
       m+3, 1];
cnc = repmat(cnc, 1, m);
cnc += kron([0:m-1], [1,1]);
cnc = repmat(cnc, 1, n);
cnc += kron([0:m+1:m*n], ones(3,m*2));

trimesh(cnc', nodes(1,:), nodes(2,:), nodes(3,:));
axis equal

midpoint = squeeze(sum(reshape(nodes(:,cnc), 3,3,[]), 2))/3;

num_elements = size(cnc, 2);
num_nodes = size(nodes, 2);

X1 = nodes(:, cnc(1, :));
X2 = nodes(:, cnc(2, :));
X3 = nodes(:, cnc(3, :));
Area = 0.5 * norm(cross(X2-X1, X3-X1, 1), 2, 'cols');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix and rhs assembly                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0_func = @(x) [0.5 * (B0y * x(3,:) - B0z * x(2,:)); ...
                0.5 * (B0z * x(1,:) - B0x * x(3,:)); ...
                0.5 * (B0x * x(2,:) - B0y * x(1,:))];
A0 = A0_func(midpoint);

Curlw = nan(3, num_elements, 3);
for ii = 1:num_elements
  for mm = 1:3
    Curlw(:, ii, mm) = curlw(mm, X1(:, ii), X2(:, ii), X3(:, ii));
  end
end

M = zeros(num_nodes);
rhs = zeros(num_nodes, 1);

fprintf('Matrix assembly\n')
tic
GG = tanh((1+1i)*thickness/2/skindepth)/((1+1i)*thickness/2/skindepth);
for ii = 1:num_elements

  [gpi, gwi] = tri_gauss_points(3, X1(:, ii), X2(:, ii), X3(:, ii));
  gpi = squeeze(gpi);
  gwi = squeeze(gwi);

  % resistance matrix
  %
  %                      /  ->         ->
  % R = 1/(GG . sigma) * | CurlT(x) . CurlT'(x) dx
  %                      /
  %                      S
  
  for mm = 1:3
    for nn = 1:3
      M(cnc(mm, ii), cnc(nn, ii)) += Curlw(:,ii,nn).' * Curlw(:,ii,mm) * Area(ii) / sigma / GG;
    end
  end

  % inductance matrix
  %
  %                                 / /   ->         ->
  % L = (mu0 * thickness)/(4*pi) *  | | (CurlT(x) . CurlT'(y)) / |x - y| dxdy
  %                                 / /
  %                                 S S
  G = zeros(numel(gwi), num_elements);
  for kk = 1:numel(gwi)
    G(kk,:) = int_Gconst_tri_mesh(nodes, cnc, gpi(:,kk));
  end
  for mm = 1:3
    for nn = 1:3
      MM = (gwi.' * G) .* ((Curlw(:,ii,mm).' * Curlw(:,:,nn))) ...
         * thickness * 1i * omega * mu0 / (4 * pi);
      M(cnc(mm, ii), :) += accumarray(cnc(nn,:).', MM, [1,num_nodes]);
    end
  end % M = R + 1i * omega * L

  % right-hand side
  %
  %                     /  ->       ->
  % rhs = -1i * omega * | A0(x) . CurlT(x) dx
  %                     /
  %                     S
  for mm = 1:3
    rhs(cnc(mm, ii)) += -1i * omega * Curlw(:,ii,mm).' * A0(:, ii) * Area(ii);
  end

end
toc
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% homogeneous dirichlet boundary conditions on T                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_left  = [1: (m+1): (m+1)*(n+1)];
idx_right = [m+1: (m+1): (m+1)*(n+1)];
idx_lower = [1: 1: (m+1)];
idx_upper = [(m+1)*n+1: 1: (m+1)*(n+1)];
idx_boundary = [idx_lower, idx_upper, idx_left, idx_right];
M(idx_boundary,:) = 0;
for ii = idx_boundary
  M(ii, ii) = 1;
end
rhs(idx_boundary) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear system solution                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diagonal preconditioner
Pinv = 1./ diag(M);
M = Pinv .* M;
rhs = Pinv .* rhs;

T = M \ rhs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing and plots                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = zeros(3, num_elements);
for ii = 1:num_elements
  for mm = 1:3
    J(:,ii) += T(cnc(mm, ii))*Curlw(:,ii,mm);
  end
end

np = 4;


Btot = zeros(3,num_elements);

prec_in_small_tria = 3;
prec_in_near_tria = 4;
prec_in_far_tria = 1;
distance_near = 0.004;
disp("start computation of Btotal\ncomputed elements:")
elem_index_vector = 1:num_elements;
pkg load parallel
par_fun = @(id_pr) parallel_losses(id_pr,np,prec_in_small_tria,...
            distance_near,prec_in_near_tria,prec_in_far_tria,...
            midpoint,num_elements,J,X1,X2,X3);

%{
  
for ii = 1:num_elements
  if mod(ii,50) == 0
    disp(ii)
  endif
  for jj = 1:num_elements
    if jj == ii
      x1tr = X1(:,ii);
      x2tr = x1tr + (X2(:,ii) - x1tr)/3;
      x3tr = x1tr + (X3(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
      x1tr = X2(:,ii);
      x2tr = x1tr + (X3(:,ii) - x1tr)/3;
      x3tr = x1tr + (X1(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
      x1tr = X3(:,ii);
      x2tr = x1tr + (X1(:,ii) - x1tr)/3;
      x3tr = x1tr + (X2(:,ii) - x1tr)/3;
      [gpst,gwst] = tri_gauss_points(prec_in_small_tria,x1tr,x2tr,x3tr);
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    elseif norm(midpoint(ii) - midpoint(jj)) < distance_near
      [gpst,gwst] = tri_gauss_points(prec_in_near_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    else
      [gpst,gwst] = tri_gauss_points(prec_in_far_tria,X1(:,jj),X2(:,jj),X3(:,jj));
      gps = squeeze(gpst);
      gws = squeeze(gwst);
      for kk = 1:size(gps,2)
        Btot(:,ii) += gws(kk).*cross(J(:,ii),gps(:,kk)-midpoint(:,ii))./(norm(gps(:,kk)-midpoint(:,ii))^3);
      end
    endif
  end
end
%}
[Btot1,Btot2,Btot3,Btot3] = parcellfun(np,par_fun,{0,1,2,3})
Btot = [Btot1,Btot2,Btot3,Btot3];
Btot = Btot.*(GG*mu0/4/pi);
Btot += [B0x,B0y,B0z]'.*ones(3,num_elements);
Bm = Btot - (dot(Btot,cross(X2-X1,X3-X1))./((norm(cross(X2-X1,X3-X1),'cols')).^2)) .* cross(X2-X1,X3-X1);


losses = thickness*real(dot(J,J)/GG/sigma - 1i*omega*dot(Bm,Bm)/GG/mu0);
losses_tot = dot(losses,Area)




%losses = 0.5 * dot(J,J,1) / sigma / GG;
%losses_tot = losses * Area.' * thickness
%fprintf('Joule losses = %e W\n', losses_tot)

Jabs_real = norm(real(J), 2, 'cols');
Jabs_imag = norm(imag(J), 2, 'cols');


figure
hold on
axis equal
hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), Jabs_real);
set(hh,'FaceColor','flat');
colormap('jet')
colorbar
quiver3(midpoint(1,:), midpoint(2,:), midpoint(3,:), ...
        real(J(1,:)), real(J(2,:)), real(J(3,:)), ...
        'k', 'linewidth', 1.5)
title('Volume current density - real part [A/m^2]')

figure
hold on
axis equal
hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), Jabs_imag);
set(hh,'FaceColor','flat');
colormap('jet')
colorbar
quiver3(midpoint(1,:), midpoint(2,:), midpoint(3,:), ...
        imag(J(1,:)), imag(J(2,:)), imag(J(3,:)), ...
        'k')
title('Volume current density - imaginary part [A/m^2]')

figure
hold on
axis equal
hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), losses);
set(hh,'FaceColor','flat');
colormap('jet')
colorbar
title('Volumetric loss density [W/m^3]')
