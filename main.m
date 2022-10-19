% clearvars
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 4.6729e6;
freq = 50;
thickness = 3e-3;
omega = 2 * pi * freq;
mu0 = 4 * pi * 1e-7;

% uniform background magnetic field (peak values)
B0x = 0;
B0y = 0;
B0z = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh construction                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0.1;
b = 0.1;
m = 24;
n = 24;

[Y,X] = meshgrid(linspace(-b/2,b/2,n+1), linspace(-a/2,a/2,m+1));
[Y,X] = meshgrid(b/2*sin(linspace(-pi/2,pi/2,n+1)),a/2*sin(linspace(-pi/2,pi/2,m+1)));

Z = zeros(size(X));
Z = 0 * (X.^2 - Y.^2);


% plot_options = 2; %1 to plot separately real and imag, 2 to plot together

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

if plot_options
  figure
  trimesh(cnc', nodes(1,:), nodes(2,:), nodes(3,:));
  axis equal
endif

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
for ii = 1:num_elements

  [gpi, gwi] = tri_gauss_points(3, X1(:, ii), X2(:, ii), X3(:, ii));
  gpi = squeeze(gpi);
  gwi = squeeze(gwi);

  % resistance matrix
  %
  %               /  ->         ->
  % R = 1/sigma * | CurlT(x) . CurlT'(x) dx
  %               /
  %               S
  for mm = 1:3
    for nn = 1:3
      M(cnc(mm, ii), cnc(nn, ii)) += Curlw(:,ii,nn).' * Curlw(:,ii,mm) * Area(ii) / sigma;
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

losses = 0.5 * dot(J,J,1) / sigma;
losses_tot = losses * Area.' * thickness;
fprintf('Joule losses = %e W\n', losses_tot)

Jabs_real = norm(real(J), 2, 'cols');
Jabs_imag = norm(imag(J), 2, 'cols');


if plot_options
  if plot_options == 1
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), Jabs_real);
    set(hh,'FaceColor','flat');
    colormap('jet')
    colorbar
    quiver3(midpoint(1,:), midpoint(2,:), midpoint(3,:), ...
            real(J(1,:)), real(J(2,:)), real(J(3,:)), ...
            'k', 'linewidth', 1.5);
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
            'k');
    title('Volume current density - imaginary part [A/m^2]')
  endif
  if plot_options == 2
    figure
    hold on
    axis equal
    hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), norm([Jabs_real;Jabs_imag],'cols'));
    set(hh,'FaceColor','flat');
    colormap('jet')
    colorbar

    title('Volume current density [A/m^2]')
  endif
  
  figure
  hold on
  axis equal
  hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), losses);
  set(hh,'FaceColor','flat');
  colormap('jet')
  colorbar
  title('Volumetric loss density [W/m^3]')
endif