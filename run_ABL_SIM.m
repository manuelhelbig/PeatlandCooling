%%% run ABL simulations to assess peatland impacts on regional climate %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig
% contact: mhelbig85@gmail.com
% edited: 2020-04-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load albedo data for peatland and forest sites
% load information on ecosystem type of sites

% define constants
ep = 0.98;          % emissivity         
sigma = 5.67E-08 ;  % Stefan Boltzmann constant      
Cp = 1005;          % Specific heat of air, J kg-1 K-1        
Rstar = 8.3144;     % universal gas constant, J mol-1 K-1
Mair = 29;          % molecular mass of air, g mol-1   
dLdT = -2370;       % derivative of latent heat of evaporation with temperature
del_t=3600;         % time step, s

% load diurnal ground heat flux and aerdoynamic conductance
% ecosystem type information [PTL == 1 (peatland), PTL == 2 (forest)]

%% load balloon sounding information (initial conditions and lapse rates)
% TEMP - initial air temperature (degC)
% MIXING_R - initial water vapor mixing ration (g kg-1)
% remove unrealistic initial value
TEMP(abs(TEMP)>40 | TEMP<0)=NaN;
MIXING_R(MIXING_R<3)=NaN;
%% run atmospheric boundary layer model
% choose ecosystem type: peatland = 1; forest = 2;
ECO=1;
% LAT - latitude of sites

% create matrices for model output
Tp_Max=NaN(300,8);      % mixed-layer potential temp (K) [4pm]
q_Max=NaN(300,8);       % mixed-layer specific humidity (g kg-1) [4pm]
h_Max=NaN(300,8);       % mixed-layer ABL height (m) [4pm]
h_Max_tot=NaN(300,8);   % maximum ABL height
Bow=NaN(300,8);         % Bowen ratio [4 pm]
VPD_all=NaN(300,8);     % vapor pressure deficit (kPa) [4pm]
esat_m=NaN(300,8);      % saturation vapor pressure
ea_m=NaN(300,8);        % actual vapor pressure

% simulate diurnal air temperature and humidity development for eight sites
% j_l - number of sites
% par_FREE_ATMOS - temp lapse rate in free atmosphere
% par_H2O_FREE_ATMOS - specific humidity lapse rate in free atmosphere
% PTL_l - number of peatland sites
% FOR_l - number of forest sites
% ALB_PTL - monthly albedo matrix for all peatland sites
% ALB_FOR - monthly albedo matrix for all forest sites
% G_FR_CL - diurnal variation in ground heat flux (forest)
% G_PT_CL - diurnal variation in ground heat flux (peatland)
% Ga_FR_CL - diurnal variation in aerodynamic conductance (forest)
% Ga_PT_CL - diurnal variation in aerodynamic conductance (peatland)

for j=1:j_l;
    j
    lat=LAT(j);
    % get lapse rates (derived from soundings)
    % only use data when lapse rates are well defined (R2 > 0.9)
    par_FREE_ATMOS{j}(1,find(R2(j,:)<0.9))=NaN;
    par_H2O_FREE_ATMOS{j}(1,find(R2(j,:)<0.9))=NaN;
    
    % only use data if valid data for all initial conditions is available
    ind=find(~isnan(MIXING_R(j,:)) & ~isnan(TEMP(j,:))...
        & ~isnan(par_FREE_ATMOS{j}(1,:)) & ~isnan(par_H2O_FREE_ATMOS{j}(1,:)));
    
    % lapse rates for potential temperature
    gammatheta_ALL=par_FREE_ATMOS{j}(1,ind);
    % lapse rates for vapor mixing ratio
    gammaq_ALL=par_H2O_FREE_ATMOS{j}(1,ind)./1000; % convert to kg kg-1
    
    % get initial air temperature and vapor mixing ratio
    TEMP_INI=TEMP(j,ind);
    MR_INI=MIXING_R(j,ind);
    
    % define parameter
    h          = 50;         % initial ABL height [m]
    dtheta     = 1.5;        % initial temperature jump at h [K]
    dq         = -0.0015;    % initial specific humidity jump at h [kg kg-1] 
    dz_h       = 150;        % Transition layer thickness [m]
    ustar      = 0.3;        % surface friction velocity [m s-1]
    
    % day of year of simulation
    DOY = 197;
    % simulate all days with available data
    for s=1:length(gammaq_ALL);
        
        % July albedo (peatland)
        if ECO==1 & s==1;
            % choose randomly one site for albedo values
            i=randsample(PTL_l,length(gammaq_ALL),true);
        % July albedo (forest) 
        elseif ECO==2 & s==1;
            i=randsample(FOR_l,length(gammaq_ALL),true);
        end
        
        if ECO==1
            met.albedo = ALB_PTL(7,i(s));
        elseif ECO==2
            met.albedo = ALB_ENF(7,i(s));
        end
        
        % initial moisture and temperature
        Q_S=MR_INI(s)./1000;                % vapor mixing ratio
        met.Tair_K = TEMP_INI(s)+273.15;    % convert to K
        theta_S = met.Tair_K; 
        
        % lapse rates (temp and moisture)
        gammatheta=gammatheta_ALL(s);
        gammaq=gammaq_ALL(s);

        % initial VPD [Pa]
        met.P_Pa = 101.3.*1000;
        met.eair_Pa = Q_S * met.P_Pa / (0.378 * Q_S + 0.622);
        esat = 100.*exp(54.8781919 - 6790.4985./ met.Tair_K...
            - 5.02808.* log(met.Tair_K));
        VPD = (esat-met.eair_Pa)./1000;

        % aerodynamic conductance (clear-sky) & ground heat flux
        Ga_CL_MIN=NaN(24*60-1,1);
        if ECO==1
            % median diurnal cycle across all peatland sites
            GA=nanmedian(Ga_PT_CL');
            energy.Gsoil = nanmedian(G_PT_CL');
        elseif ECO==2;
            % median diurnal cycle across all forest sites
            GA=nanmedian(Ga_FR_CL');
            energy.Gsoil = nanmedian(G_FR_CL');
        end
        
        % create continuous 1-min time series for model simulations
        for k=1:24*60-1;
            if k<23*60
                Ga_CL_MIN(k)=GA(floor(k/60)+1)...
                    +(GA(floor(k/60)+2)-GA(floor(k/60)+1))/60.*(k-(floor(k/60)*60));
            elseif k == 23*60;
                Ga_CL_MIN(k)=GA(floor(k/60));
            else
                Ga_CL_MIN(k)=GA(floor(k/60)+1)...
                    +(GA(1)-GA(floor(k/60)+1))/60.*(k-(floor(k/60)*60));
            end
        end
        GA=Ga_CL_MIN;
        
        % run ABL model
        theta_ts=[];    % pot temp
        h_ts=[];        % ABL height
        ts_ts=[];       % surf temp
        q_ts=[];        % mixing ratio
        H_ts=[];        % sensible heat flux
        LE_ts=[];       % latent heat flux
        LE_PM_ts=[];    % latent heat flux (from Penman-Monteith)
        VPD_ts=[];      % vapor pressure deficit
        AV_ts=[];       % available energy flux
        
        % run for 14h (05h-19h)
        % half-hour time steps
        for v=1:14*2;
            % hour of day
            HH=(v-1)./2+5;
        
            if v == 1
                % run ABL model (adjusted from CLASS MXL model) for first half hour
                [h_out,theta_out,dtheta_out,q_out,dq_out,dz_h_out,H_out,LE_out,LE_PM_out,AV_out,ts_out]...
                    = ABL_MODEL_v1(h,dtheta,gammatheta,dq,gammaq,dz_h,ustar,theta_S,Q_S,...
                    met,energy,GA,nanmedian(par_gs(PTL==ECO,:)),...
                    nanmedian(par_gs_SWIN(PTL==ECO,:)),nanmedian(GS_max(PTL==ECO)),...
                    HH,DOY+HH./24,lat);
                theta_S_test(s) = theta_S;
                gammatheta_test(s) = gammatheta;
            else
                 % run ABL model (adjusted from CLASS MXL model) for half hour (use last model output to run model)
                [h_out,theta_out,dtheta_out,q_out,dq_out,dz_h_out,H_out,LE_out,LE_PM_out,AV_out,ts_out]...
                    = ABL_MODEL_v1(h_out(end),dtheta_out(end),gammatheta,dq_out(end),...
                    gammaq,dz_h_out(end),ustar,theta_out(end),q_out(end),...
                    met,energy,GA,nanmedian(par_gs(PTL==ECO,:)),...
                    nanmedian(par_gs_SWIN(PTL==ECO,:)),nanmedian(GS_max(PTL==ECO)),...
                    HH,DOY+HH./24,lat);
            end
            % create half-hour time series (as mean of 1-min model output)
            theta_ts=horzcat(theta_ts,nanmean(theta_out));
            h_ts=horzcat(h_ts,nanmean(h_out));
            ts_ts=horzcat(ts_ts,nanmean(ts_out));
            q_ts=horzcat(q_ts,nanmean(q_out));
            H_ts=horzcat(H_ts,nanmean(H_out));
            LE_ts=horzcat(LE_ts,nanmean(LE_out));
            LE_PM_ts=horzcat(LE_PM_ts,nanmean(LE_PM_out));
            AV_ts=horzcat(AV_ts,nanmean(AV_out));

            eair_Pa = nanmean(q_out) * met.P_Pa / (0.378 * nanmean(q_out) + 0.622);
            esat = 100.*exp(54.8781919 - 6790.4985./ nanmean(theta_out)...
            - 5.02808.* log(nanmean(theta_out)));
            VPD_ts = horzcat(VPD_ts,(esat-eair_Pa)./1000);
            
            % output of afternoon conditions
            if HH == 16;
                if isreal(theta_ts)
                    Tp_Max(s,j) = nanmean(theta_out);
                    q_Max(s,j) = nanmean(q_out);
                    h_Max(s,j) = nanmean(h_out);
                    Bow(s,j) = nanmean(H_out)./nanmean(LE_out);
                    VPD_all(s,j) = (esat-eair_Pa)./1000;
                    esat_m(s,j) = esat./1000;
                    ea_m(s,j) = eair_Pa./1000;              
                end
            end
            
            if HH == 18.5;
                h_Max_tot(s,j) = max(h_ts);
            end
        end   
    end
        % filter for unrealistic model runs (based on difference between max temp and initial temp
        diff=real(theta_S_test'-Tp_Max(1:length(gammaq_ALL),j));
        % remove upper and lower 2.5 percent and if diurnal variation in Tpot > 25K
        index=find(diff<prctile(diff,2.5) | diff>prctile(diff,97.5) | abs(diff)>25);
        Tp_Max(index,j) = NaN;
        q_Max(index,j) = NaN;
        h_Max(index,j) = NaN;
        h_Max_tot(index,j) = NaN;
        Bow(index,j) = NaN;
        VPD_all(index,j) = NaN;
        esat_m(index,j) = NaN;
        ea_m(index,j) = NaN;
        clear theta_S_test index
end