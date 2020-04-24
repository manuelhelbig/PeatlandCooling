function [H, LE, LE_PM, TS, AVAIL]...
    = SURF_EN_BALANCE(met,energy,GA,par_VPD,par_SWIN,GS_MAX,HH,DOY,lat)
%%%%%%%%%%%%%%%%%%%% surface energy balance model %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig
% McMaster University
% last edit: 2020-04-23
% email: mhelbig85@gmail.com
% adapted code from Dennis Baldocchi's C+ code (UC Berkeley)
% https://nature.berkeley.edu/biometlab/BiometWeb/ET_pbl.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
    % met.Tair_K = near-surface air temperature [K]
    % met.albedo = albedo
    % met.eair_Pa = vapor pressure [Pa]
    % met.RgInMax = maximum irradiance [W m-2]
    % met.P_Pa = barometric pressure [Pa]
    % energy.Gsoil = ground heat flux
    % GA          aerodynamic conductance [m s-1]
    % par_VPD     multiple constraint function of gs for VPD
    % par_SWIN    multiple constraint function of gs for light
    % GS_MAX      maximum surface conductance
    % HH          hour of day
    % DOY         day of year
    % lat         latitude

% outputs
    % H           sensible heat flux [W m-2]
    % LE          latent heat flux [W m-2]
    % LE_PM       latent heat flux from Penman-Monteith [W m-2]
    % TS          surface temperature
    % AVAIL       available energy flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ep = 0.98;          % emissivity (perfect radiatior)

% clear-sky emissivity (Brutsaert 1975)
A=1.24;
B=1/7;
ep_cs = A*((met.eair_Pa/100)./met.Tair_K)^B;

% define constants
sigma = 5.67E-08 ;  % Stefan Boltzmann constant      
Cp = 1005;          % Specific heat of air, J kg-1 K-1        
Rstar = 8.3144;     % universal gas constant, J mol-1 K-1
Mair = 29;          % molecular mass of air, g mol-1   
dLdT = -2370;       % derivative of latent heat of evaporation with temperature

% convert from Pa to kPa
met.P_kPa = met.P_Pa./1000;

% aerodnamic conductance
Ga_INI=GA;

% calculate solar radiation based on latitude [W m-2] (Allen et al 2005)
% 5am to 7pm
S0 = 1368;
t_day=(DOY-ceil(DOY)).*24;
DOY=ceil(DOY);
sda    = 0.409 * cos(2.*pi()* (DOY - 173)/365);
sinlea = sin(2. * pi()* lat/360) * sin(sda) - cos(2.* pi() * lat./ 360)*...
        cos(sda) * cos(2.* pi() * (t_day * 3600)./86400); % - 2.* pi() * lon/ 360);
sinlea = max(sinlea, 0.0001);
met.RgIn  = (S0 * sinlea).*0.8;
clear sinlea sda t_day DOY S0

% Infrared radiation from sky, W m-2, using algorithm from Norman [W m-2]
met.LongIn = ep_cs * sigma * met.Tair_K.^4; 
energy.Lout = ep_cs * sigma * met.Tair_K.^4;

% incoming radiation (SW+LW-G)
met.Q_In=met.RgIn.*(1-met.albedo)+ep.*met.LongIn-energy.Gsoil;

% saturation vapor pressure [Pa]
esat = 100.*exp(54.8781919 - 6790.4985./ met.Tair_K...
    - 5.02808.* log(met.Tair_K));
% vapor pressure deficit [Pa]
met.vpd_Pa=esat-met.eair_Pa;
% convert Pa to kPa
VPD=met.vpd_Pa./1000;

%%%%%%%%%%%%% Resistances and Conductances %%%%%%%%%%%%%%%%%%%
SWIN = met.RgIn;
% VPD modifier for surface conductance [range: 0-1]
gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));
fVPD=gsFxn(par_VPD,VPD);
% SWin modifier  [range: 0-1]
gsFxnSWIN=@(params,SWIN) ((params(1).*params(2).*SWIN)./(params(1)+params(2).*SWIN));
fSWin=gsFxnSWIN(par_SWIN,SWIN);

% if modifier is 0
if fSWin==0;
    fSWin=0.01;
end

if fVPD==0;
    fVPD=0.01;
end

% calculate surface conductance
var.Gs=GS_MAX.*fVPD.*fSWin;
var.Rs=1./var.Gs; % surface resistance
var.Rheat=1./Ga_INI; % aerodynamic resistance
var.Gh=1./var.Rheat;
var.Gw=1./(var.Rheat+var.Rs);

% LATENT HEAT OF VAPORIZATION, J KG-1
lmbda = 3149000 - 2370 * met.Tair_K;
var.psych = met.P_Pa * Cp / (0.622 * lmbda);    % Pa K-1
% Slope of the saturation vapor pressure-temperature curve
% Thermodynamic relationship from HESS    (Pa K-1)
esat=100.*exp(54.8781919 - 6790.4985./ met.Tair_K - 5.02808.* log(met.Tair_K));
var.slope = esat * lmbda * 18./ (Rstar * met.Tair_K * met.Tair_K * 1000);

% Initialize Derived Quantities

% for convenience, Gsoil is substracted from Qin
met.Q_In=met.RgIn.*(1-met.albedo)+ep*met.LongIn;

energy.Rnet= met.RgIn*(1-met.albedo)+ep*met.LongIn-energy.Lout;
energy.Avail = energy.Rnet - energy.Gsoil(ceil(HH));

var.air_density = met.P_kPa * Mair / (Rstar * met.Tair_K); 
energy.Rniso=met.RgIn*(1.-met.albedo)+ep*met.LongIn-ep*sigma*met.Tair_K.^4;
var.Gr=4.*ep*sigma*met.Tair_K.^3/(var.air_density*Cp);

% energy balance computations using quadradic solutions
% to the surface energy balance (both for Tsfc and LE)
    
% Quadratic solution of LE at canopy scale

% aLE^2 + bLE + c = 0
    
dest = esat * lmbda * 18. / (Rstar * met.Tair_K * met.Tair_K * 1000);
d2est = -2 * esat * lmbda * 18./ (Rstar * met.Tair_K * met.Tair_K * met.Tair_K * 1000)...
    +  dest * lmbda * 18./ (Rstar * met.Tair_K * met.Tair_K * 1000);

energy.Lout_Ta = ep_cs * sigma * met.Tair_K.^4;
lecoef=var.air_density*0.622*lmbda*var.Gw/(met.P_Pa);
hcoef=var.air_density*Cp*var.Gh;
longcoef=4*ep*sigma*met.Tair_K.^3;

Acoef=lecoef*d2est/(2*(hcoef+longcoef));

Bcoef=-(hcoef+longcoef)-lecoef*dest + Acoef*(2.*energy.Lout_Ta +2*energy.Gsoil(ceil(HH))- 2.*met.Q_In);

Ccoef=(hcoef+longcoef)*lecoef*met.vpd_Pa+lecoef*dest*(met.Q_In - energy.Lout_Ta -energy.Gsoil(ceil(HH))) + ...
			Acoef*(met.Q_In*met.Q_In+energy.Lout*energy.Lout_Ta + energy.Gsoil(ceil(HH))*energy.Gsoil(ceil(HH))-...
			2.*energy.Lout_Ta*met.Q_In-2*met.Q_In*energy.Gsoil(ceil(HH))+2*energy.Lout_Ta*energy.Gsoil(ceil(HH)));
% solve for LE
% a LE^2 + bLE + c = 0

le1 = (-Bcoef + (Bcoef *Bcoef - 4.* Acoef * Ccoef).^0.5) / (2. * Acoef);
le2 = (-Bcoef - (Bcoef * Bcoef - 4. * Acoef * Ccoef).^0.5) / (2. * Acoef);

energy.LE=le2;

%solve for Ts
% aT^2 + bT + c =0
aT = 6.* ep * sigma * met.Tair_K * met.Tair_K + d2est * lecoef/2;
bT = longcoef + hcoef + dest * lecoef;
cT =  energy.Lout_Ta + lecoef * met.vpd_Pa - met.Q_In + energy.Gsoil(ceil(HH));

del_Tk=(-bT + (bT * bT - 4. * aT * cT).^0.5) / (2. * aT);
met.Tsfc_K= met.Tair_K + del_Tk;

% H IS SENSIBLE HEAT FLUX
energy.H = hcoef * (del_Tk);
H=energy.H;
LE=energy.LE;

%LOUT IS LONGWAVE EMITTED RADIATION
energy.Lout = ep * sigma * met.Tsfc_K.^4;

% Compute isothermal Penman-Monteith LE and equilibrium LE
    
energy.Rnet= met.RgIn*(1-met.albedo)+ep*met.LongIn-energy.Lout;
energy.Avail = energy.Rnet - energy.Gsoil(ceil(HH));
  
energy.LE_iso = (var.slope * (energy.Rniso-energy.Gsoil(ceil(HH))) + Cp * var.air_density * met.vpd_Pa*(var.Gr+var.Gh))/...
(var.slope + var.psych * (var.Gr+var.Gh)/var.Gw);
    
energy.LE_PM = (var.slope * (energy.Avail) + Cp * var.air_density * met.vpd_Pa*var.Gh)/...
    (var.slope + var.psych * var.Gh/var.Gw);

energy.LE_eq = var.slope * energy.Avail / (var.slope + var.psych);

LE_PM = energy.LE_PM;
AVAIL = energy.Avail;
TS = met.Tsfc_K;