%% alpha-Fe
% parameters for bulk diffusion 
pre_exponential = 2e-4; % m^3/s
active_energy = 251*1e3; % kJ/mol  2.6015eV
R = 8.314; % J/mol/K
Boltz    = 1.3807e-23; % J/K
T=750; % K
Dv = pre_exponential*exp(-active_energy/R/T); %

% parameters for core diffusion 
pre_exponential = 1e-23; % m^4/s
active_energy = 174*1e3; % kJ/mol  1.8034eV
R = 8.314; % J/mol/K
Boltz    = 1.3807e-23; % J/K
T=750; % K
D_core = pre_exponential*exp(-active_energy/R/T);  % core diffusivity
Atomvol  = (2.856)^3*1e-30;  %%%%%%%       volum of the atom
Dc=D_core*Atomvol/Boltz/T;


%% parameters from scientific report Thomas, 2016
Boltz    = 1.3807e-23; % J/K
T=750; % K
Dsc=1.0e-19; % m^2/s
Atomvol  = (2.856)^3*1e-30;  % volum of the atom
D_core = Dsc*pi*(0.247e-9)^2*(3.8e-9)^3;  % ac×Dc*V

Dc=D_core/Boltz/T;


% parameters for climb

Coeff_core0 = 1.18e-5;
formeng  = 0.67*1.6022e-19;%%%%%   formation energy of vacancy
migeng   =0.6* 0.68*1.6022e-19; %        migration energy of vacancy
Boltz    = 1.3807e-23;
Coeff_core =Coeff_core0*exp(-migeng/Boltz/tem);  % core diffusivity;

Atomvol  = 16.3e-30;  %%%%%%%       volum of the atom
rc       = 5*bur*1e-9;   % radius of the dislocation core area

Dc = 1e-3*Coeff_core*(2*pi*rc^2)*Atomvol/Boltz/tem/(bur*1e-9)^5;