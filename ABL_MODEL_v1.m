function [h_out,theta_out,dtheta_out,q_out,dq_out,dz_h_out,H_out,LE_out,LE_PM_out,AE_out,Ts_out]...
    = ABL_MODEL_v1(h,dtheta,gammatheta,dq,gammaq,dz_h,ustar,theta_S,Q_S,...
    met,energy,GA,par_VPD,par_SWin,GS_MAX,HH,DOY,lat)
%%%%%%%%%%%%%%%%%%%%% ABL slab model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig
% McMaster University
% email: mhelbig85@gmail.com
% adapted code from CLASS MXL model
% see for more information: https://classmodel.github.io/
% CLASS model was developed by the Meteorology and Air Quality section, Wageningen University and Research
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
    % h           initial ABL height
    % dtheta      initial temperature jump at h [K]
    % gammatheta  free atmosphere potential temperature lapse rate [K m-1]
    % dq          initial specific humidity jump at h [kg kg-1]
    % gammaq      free atmosphere specific humidity lapse rate [kg kg-1 m-1]
    % dz_h        transition layer thickness [m]
    % ustar       friction velocity [m s-1]
    % theta_S     initial mixed-layer potential temperature [K]
    % Q_S         initial mixed-layer specific humidity [kg kg-1]
    % met         structure with air temp [K], atmospheric pressure [Pa]
    % energy      structure with...
    % GA          aerodynamic conductance [m s-1]
    % par_VPD     multiple constraint function of gs for VPD
    % par_SWIN    multiple constraint function of gs for light
    % GS_MAX      maximum surface conductance
    % HH          hour of day
    % DOY         day of year
    % lat         latitude

% outputs
    % h_out
    % theta_out
    % dtheta_out
    %q_out
    % dq_out 
    % dz_h_out
    % H_out
    % LE_out
    % LE_PM_out
    % AE_out
    % Ts_out

% initial mixed-layer potential temperature [K]
met.Tair_K = theta_S;

% atmospheric pressure [Pa]
met.P_Pa = 101.3.*1000;

% partial pressure of water vapor [Pa]
met.eair_Pa = Q_S * met.P_Pa / (0.378 * Q_S + 0.622);

% aerodynamic conductance [m s-1]
GA_INI = GA(HH*60);

% run surface energy balance scheme to derive H and LE
[H, LE, LE_PM, TS, AVAIL]= SURF_EN_BALANCE(met,energy,GA_INI,par_VPD,par_SWin,...
    GS_MAX,HH,DOY,lat);

% initialize constants
Lv         = 2.5e6;                 % heat of vaporization [J kg-1]
cp         = 1005;                 % specific heat of dry air [J kg-1 K-1]
rho        = 1.2;                   % density of air [kg m-3]
k          = 0.4;                   % Von Karman constant [-]
g          = 9.81;                  % gravity acceleration [m s-2]
Rd         = 287;                  % gas constant for dry air [J kg-1 K-1]
Rv         = 461.5;                 % gas constant for moist air [J kg-1 K-1]
bolz       = 5.67e-8;               % Bolzman constant [-]
rhow       = 1000;                 % density of water [kg m-3]
S0         = 1368;                 % solar constant [W m-2]

% initialize mixed layer potential temp
theta=theta_S(1);
% advection of heat [K s-1]
advtheta   = 0;       
% entrainment ratio for virtual heat [-]
beta       = 0.2;               

% initial sensible heat flux
wtheta   = H(1)  / (rho * cp);  % surface kinematic heat flux [K m s-1]   
wthetae    = [];                % entrainment kinematic heat flux [K m s-1]

wstar      = 0;                 % convective velocity scale [m s-1]

% moisture 
q          = Q_S(1);   % initial mixed-layer specific humidity [kg kg-1]
advq       = 0;       % advection of moisture [kg kg-1 s-1]

% initial latent heat flux
wq       = LE(1) / (rho * Lv);   % surface kinematic moisture flux [kg kg-1 m s-1]
wqe        = 0;                  % entrainment moisture flux [kg kg-1 m s-1]
%wqM        = 0;                  % moisture cumulus mass flux [kg kg-1 m s-1]

% mixed-layer input
Ps         = 101300;   % surface pressure [Pa]
divU       = 0;        % horizontal large-scale divergence of wind [s-1]
fc         = 1.e-4;    % Coriolis parameter [m s-1]

% initialize cumulus parameterization
sw_cu      = 0;         % Cumulus parameterization switch
sw_wind    = 0;         % prognostic wind switch
ac         = 0;         % Cloud core fraction [-]
M          = 0;         % Cloud core mass flux [m s-1] 
wqM        = 0;         % Cloud core moisture flux [kg kg-1 m s-1] 
wqe        = 0;         % entrainment moisture flux [kg kg-1 m s-1]

u          = 6;         % initial mixed-layer u-wind speed [m s-1]
du         = 4;         % initial u-wind jump at h [m s-1]
gammau     = 0;         % free atmosphere u-wind speed lapse rate [s-1]
advu       = 0;         % advection of u-wind [m s-2]
v          = -4.0;      % initial mixed-layer v-wind speed [m s-1]
dv         = 4.0;       % initial u-wind jump at h [m s-1]
gammav     = 0;         % free atmosphere v-wind speed lapse rate [s-1]
advv       = 0;         % advection of v-wind [m s-2]

%z0h        = z0m./10;   % roughness length for scalars [m]

sw_fixft   = 0;         % Fix the free-troposphere switch
dFz        = 0;         % cloud top radiative divergence [W m-2]
sw_shearwe = 0;     	% shear growth mixed-layer switch

dt         = 60;       % time step [s]
runtime    = 30*60;    % run for 30 mins [s]

% initialize time variables
tsteps = floor(runtime / dt);
t      = 0;

% initial potential temp of mixed layer
TA_INI = theta_S(1)-273.15;
    
for k = 1:tsteps;

    % calculate virtual temperature and virtual kinematic heat flux
    thetav   = theta  + 0.61.*theta.*q;
    wthetav  = wtheta + 0.61.*theta.*wq;
    dthetav  = (theta + dtheta) * (1 + 0.61.* (q + dq)) - theta.* (1 + 0.61.*q);

    % Mixed-layer top properties
    P_h    = Ps - rho * g * h; % pressure at ABL top
    T_h    = theta - g/cp * h; % potential temperature at top of ABL
    
    % calculate relative humidity
    esat = 0.611e3 * exp(17.2694 * (T_h - 273.16)...
        / (T_h - 35.86));   % saturation vapor pressure
    qsat=0.622 * esat./P_h;              % saturation specific humidity
    RH_h   = q / qsat;

    % Find lifting condensation level iteratively
    if t == 0;
                lcl = h;
                RHlcl = 0.5;
    else
                RHlcl = 0.9998; 
    end

    itmax = 30;
    it = 0;
    while ((RHlcl <= 0.9999) | (RHlcl >= 1.0001)) & it<itmax;
                lcl    = lcl+ (1-RHlcl).*1000;
                p_lcl        = Ps - rho * g * lcl;
                T_lcl        = theta - g/cp * lcl;
                esat_lcl = 0.611e3 * exp(17.2694 * (T_lcl - 273.16) / (T_lcl - 35.86));
                qsat_lcl=0.622 * esat_lcl./p_lcl;
                RHlcl        = q / qsat_lcl;
                it          = it + 1;
    end

    % calculate convective velocity scale w* 
    if wthetav > 0;
        wstar = ((g * h * wthetav) / thetav)^(1/3);
    else
        wstar  = 1e-6;
    end
    
    % Virtual heat entrainment flux (as a fraction of surface flux [20%])
    wthetave    = -beta * wthetav;
    
    % compute mixed-layer tendencies
    if sw_shearwe == 1;
        we    = (-wthetave + 5 * ustar^3 * thetav / (g * h)) / dthetav;
    else
        we    = -wthetave / dthetav;
    end

    
    % Calculate mixed-layer top relative humidity variance (Neggers et. al 2006/7)
    if wthetav > 0;
        q2_h   = -(wqe  + wqM )*dq*h/(dz_h * wstar);
    else
        q2_h   = 0;
    end

    % calculate cloud core fraction (ac), mass flux (M) and moisture flux (wqM)
    esat_h = 0.611e3 * exp(17.2694 * (T_h - 273.16) / (T_h - 35.86));
    qsat_h=0.622 * esat_h./P_h;
    ac     = max(0, 0.5 + (0.36 * atan(1.55 * ((q - qsat_h) / q2_h^0.5))));
    M      = ac * wstar;
    wqM    = M * q2_h^0.5;
    
    %%%%%%%%% run mixed layer %%%%%%%%%
    % decompose ustar along the wind components
    uw = -sign(u) * (ustar^4. / (v^2. / u^2 + 1))^(0.5);
    vw = -sign(v) * (ustar^4. / (u^2. / v ^2 + 1))^(0.5);

    % calculate large-scale vertical velocity (subsidence)
    ws = -divU * h;

    % calculate compensation to fix the free troposphere in case of subsidence 
    if sw_fixft==1;
        w_th_ft  = gammatheta * ws;
        w_q_ft   = gammaq     * ws;
    else
        w_th_ft  = 0;
        w_q_ft   = 0;
    end

    % calculate mixed-layer growth due to cloud top radiative divergence
    wf = dFz / (rho * cp * dtheta);

    % Don't allow boundary layer shrinking if wtheta < 0
    if we < 0;
        we = 0;
    end
    
    % Calculate entrainment fluxes
    wthetae     = -we * dtheta; % temperature
    wqe         = -we * dq; % moisture
    htend       = we + ws + wf - M;
    % temperature & moisture at end of timestep
    thetatend   = (wtheta - wthetae) / h + advtheta;
    qtend       = (wq     - wqe     - wqM  ) / h + advq;
    
    % temperature & moisture jump at h
    dthetatend  = gammatheta * (we + wf - M) - thetatend + w_th_ft;
    dqtend      = gammaq     * (we + wf - M) - qtend     + w_q_ft;

    % assume u + du = ug, so ug - u = du
    if sw_wind==1;
        utend       = -fc * dv + (uw + we * du)  / h + advu;
        vtend       =  fc * du + (vw + we * dv)  / h + advv;

        dutend      = gammau * (we + wf - M) - utend;
        dvtend      = gammav * (we + wf - M) - vtend;
    end

    % tendency of the transition layer thickness
    if ac > 0 | (lcl - h < 300)
        dztend = ((lcl - h)-dz_h) / 7200;
    else
        dztend = 0;
    end

    % set values previous time step
    h0      = h;

    theta0  = theta;
    dtheta0 = dtheta;
    q0      = q;
    dq0     = dq;

    u0      = u;
    du0     = du;
    v0      = v;
    dv0     = dv;

    dz0     = dz_h;

    % integrate mixed-layer equations
    h        = h0      + dt * htend;
    theta    = theta0  + dt * thetatend;
    dtheta   = dtheta0 + dt * dthetatend;
    q        = q0      + dt * qtend;
    dq       = dq0     + dt * dqtend;
    dz_h     = dz0     + dt * dztend;

    % Limit dz to minimal value
    dz0 = 50;
    if dz_h < dz0
        dz_h = dz0;
    end
    
    h_out(k)=h;
    theta_out(k)=theta;
    dtheta_out(k)=dtheta;
    q_out(k)=q;
    dq_out(k)=dq;
    dz_h_out(k)=dz_h;
    
    % calculate kinematic heat fluxes
    met.Tair_K = theta;
    met.eair_Pa = q * met.P_Pa / (0.378 * q + 0.622);
    GA_INI = GA(floor((HH+(k/60))*60));
    
    [H1, LE1, LE1_PM, TS1, AVAIL1]= SURF_EN_BALANCE(met,energy,GA_INI,par_VPD,par_SWin,...
        GS_MAX,HH+(k/60),DOY,lat);
    Ts_out(k) = TS1;
    H_out(k) = H1;
    LE_out(k) = LE1;
    LE_PM_out(k) = LE1_PM;
    AE_out(k) = AVAIL1;
    wtheta   = H1  / (rho * cp);
    wq       = LE1 / (rho * Lv);
    %{
    if k<=tsteps-1;
        wtheta   = H(k+1)  / (rho * cp);
        wq       = LE(k+1) / (rho * Lv);
    end
    %}
end
