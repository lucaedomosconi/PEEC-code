plot_options = [0,1,1,1,0,1];
% 1: plot 3d grid
% 2: plot average volume current density real and imag
% 3: plot average volume current density norm
% 4: plot volumetric loss density
% 5: plot B on line and/or compare with benchmark
% 6: plot B on plate -> needs compute_auto_induced_B = 1!!
%                       and change Bm definition
spire_generated_field = 0;
compute_auto_induced_B = 1;
G_form = 1;
postprocess = 0; % 1 to run only postprocess

a = 0.1;
b = 0.1;
m = 32;
n = 32;
mod_mesh = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu0 = 4 * pi * 1e-7;
mr = 1.; % must be 1!
mu = mu0*mr;
sigma = 4.6729e6;
freq = 50;
thickness = 3e-3;
% warning: the following values may be modified by n_main ! %%%%%%%%%%%%%%%%%%%%
omega = 2 * pi * freq;                                                         %
skindepth = sqrt(2/omega/mu/sigma);                                            %
GG = tanh((1+1i)*thickness/2/skindepth)/((1+1i)*thickness/2/skindepth);        %
if !G_form                                                                     %
  GG = 1.;                                                                     %
endif                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uniform background magnetic field (peak values)
B0x = 0;
B0y = 0;
B0z = 1;


% data
folder_str_prefix = "line_simulation_B_";
freq_str = ["50Hz";"1KHz";"2KHz";"5KHz"];
freq_id = 1;
fields_id = [1,3,4,6];

for freq = [50,1000,2000,5000]
n_main
postprocess_plots
end


%{

n_main
postprocess_plots

for freq = [50,1000,2000,5000]
n_main
postprocess_plots
end

for thickness = 1e-3:1e-3:5e-3
GG = tanh((1+1i)*thickness/2/skindepth)/((1+1i)*thickness/2/skindepth)
n_main
end

G_form = 0
n_main
BtotLine_old = BtotLine;
G_form = 1
n_main
postprocess_plots

freq = 50
freq_id = 1
G_form = 0
n_main
BtotLine_old = BtotLine;
G_form = 1
n_main
postprocess_plots
freq = 1000
freq_id = 2
G_form = 0
n_main
BtotLine_old = BtotLine;
G_form = 1
n_main
postprocess_plots
freq = 2000
freq_id = 3
G_form = 0
n_main
BtotLine_old = BtotLine;
G_form = 1
n_main
postprocess_plots
freq = 5000
freq_id = 4
G_form = 0
n_main
BtotLine_old = BtotLine;
G_form = 1
n_main
postprocess_plots


%}