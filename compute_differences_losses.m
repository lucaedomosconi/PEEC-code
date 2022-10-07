clearvars

plot_diff = 0;
plot_errors = 1;

if plot_diff

plot_options = 0;
temp_G = G_form;
G_form = 1;
n_main
J_int_n = J;
losses_n = losses;
Jabs_real_n = Jabs_real;
Jabs_imag_n = Jabs_imag;
G_form = 0;
n_main
J_int = J;
G_form = temp_G;

figure
hold on
axis equal
hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), norm(J_int_n,'cols')-norm(J_int,'cols'));
set(hh,'FaceColor','flat');
colormap('jet')
colorbar
title('Volume current density [A/m^2]')


figure
hold on
axis equal
hh = trisurf(cnc', nodes(1,:), nodes(2,:), nodes(3,:), losses_n - losses);
set(hh,'FaceColor','flat');
colormap('jet')
colorbar
title('Volumetric loss density [W/m^3]')

endif
if plot_errors
%       original    modified    Comsol
lossM = [2.416363e3 2.416025e3 2.4256e3; % 50Hz
         4.937592e5 4.811140e5 5.1913e5; % 1000Hz
         8.820190e5 8.437866e5 9.5829e5; % 2000Hz
         1.253207e6 1.172258e6 1.4137e6];% 4000Hz
figure
plot([50,1000,2000,4000],lossM(:,2),[50,1000,2000,4000],lossM(:,3))
axis([0,4000]);
legend('modified code',' Comsol benchmark')
residuals = (lossM(:,2)-lossM(:,3))./lossM(:,3)
endif
