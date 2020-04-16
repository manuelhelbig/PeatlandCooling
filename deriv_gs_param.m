%%% derive fit parameters for multiple constraint model of surface cond %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Manuel Helbig
% contact: mhelbig85@gmail.com
% edited: 2020-04-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define window of analysis
% number of of sites
lng = 96;
% select start and end date (e.g. July)
DOY_start=183;
DOY_end=214;
% select only daytime data
time_start = 5;
time_end = 9;

% define fitting options
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';

%% choose peatlands, undisturbed evergreen needleleaf forests, and others
% PFT_NUM = ecosystem type (PTL, ENF)
% DIST_NUM = disturbed (1) or undisturbed (0)
for k=1:lng;
    if strcmp(PFT_NUM{k},'PTL');
        PFT_NUMERIC(k)=1;
    elseif strcmp(PFT_NUM{k},'ENF') & DIST_NUM(k)==0;
        PFT_NUMERIC(k)=0;
    else
        PFT_NUMERIC(k)=2;
    end
end
%% derive multiple constraint function of surface conductance (VPD)
% input:
% DOY: day of year time series
% HH: hour of day time series
% VPD: vapour pressure deficit
% H: sensible heat flux
% LE: latent heat flux
% GS: surface conductance
% ga: aerodynamic conductance
for k=1:lng;
    k
    % select surface conductance to analyse
        % data for July
        % between 5h and 19h local time
        % minimum available energy flux = 25 Wm-2
    
    AV=LE+H; % available energy flux
    GS(GS<0)=NaN; % remove negative surface conductance
    
    % create indices for selection
    rem = find((DOY<DOY_start | DOY>=DOY_end) ...
            | (HH<time_start | HH>=time_end) | AV<25);
    
    VPD(rem)=NaN;
    GS(rem)=NaN;
    
    % derive maximum surface conductance
    GS_MAX=GS;
    GS_MAX(VPD<=0.5)=NaN;
    GS_max(k)=prctile(GS_MAX,97.5); % 97.5th percentile
    clear GS_MAX
    
    % get aerodynamic conductance
    ga(rem)=NaN;
    ga_SIT(k)=nanmedian(ga); % median site ga
    
    % VPD bins for surface conductance fitting
    PRC_VPD=0.1:0.1:3; % 0.1 kPa per bin
    i=find(~isnan(GS));
    VPD_UP = [];
    GS_UP = [];
    % derive upper boundary of gs
    for l =1:length(PRC_VPD)+1;
            if l == 1;
                ind=find(VPD<=PRC_VPD(l));     
            elseif l<=length(PRC_VPD);
                ind=find(VPD>PRC_VPD(l-1) & VPD<=PRC_VPD(l));  
            else
                ind=find(VPD>PRC_VPD(l-1));
            end
            
            if length(ind)>3 & ~isempty(i)
            	% select data for VPD bin
                VPD_sel = VPD(ind);
                GS_sel = GS(ind);
                
                % remove outliers using quartile method (GS)
                IQR = prctile(GS_sel,75)-prctile(GS_sel,25);
                outl = find(GS_sel>prctile(GS_sel,75)+1.5*IQR | GS_sel<prctile(GS_sel,25)-1.5*IQR);
                GS_sel(outl)=NaN;
                
                % calculate mean for VPD bins
                VPD_BINS(k,l)=nanmean(VPD_sel);
                GS_BINS(k,l)=nanmean(GS_sel);
                GS_BINS_SD(k,l)=nanstd(GS_sel);
                
                % select upper boundary data
                x1=VPD_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                y1=GS_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
              
                clear GS_sel VPD_sel RN_sel
                % calculate median for upper boundary per bin
                VPD_BINS_UP(k,l)=nanmedian(x1);
                GS_BINS_UP(k,l)=nanmedian(y1);
                
                VPD_UP=vertcat(VPD_UP,x1);
                GS_UP=vertcat(GS_UP,y1);
                clear x1 y1 v1 w1
            else
                VPD_BINS_UP(k,l)=NaN;
                GS_BINS_UP(k,l)=NaN;
            end
    end 
    % remove gs data for VPD < 0.5 kPa
    GS_BINS_UP(k,VPD_BINS_UP(k,:)<=0.5)=NaN;
    GS_UP(VPD_UP<=0.5)=NaN;
    
    % fit model to VPD-gs relationship
    
            if sum(~isnan(GS_UP))>10;
                % only fit model to VPD > 0.5kPa
                thr=0.5;
                GS_UP(VPD_UP<thr)=NaN;
                VPD_UP(VPD_UP<thr)=NaN;
                x1 = VPD_UP;
                y = (GS_UP./GS_max(k)); % normalize by gsmax
                i1 = find(~isnan(x1) & ~isnan(y));
                gsFxn=@(params,VPD) (params(2)+(1+params(1)./(sqrt(VPD))));
                
                par0_gs=[10 -5];
                par_gs_VPD(k,:) = nlinfit(x1(i1),y(i1),gsFxn,par0_gs,opts);
                
                % get goodness of fit metrics
                y=GS_BINS_UP(k,:)./GS_max(k);
                x=VPD_BINS_UP(k,:);
                RMSE(k)=sqrt(mean(y(~isnan(x) & ~isnan(y) & x>1)-(gsFxn(par_gs_VPD(k,:),x(~isnan(x) &~isnan(y) & x>1)))).^2);
                R=corrcoef(y(~isnan(x) & ~isnan(y) & x>1),(gsFxn(par_gs_VPD(k,:),x(~isnan(x) &~isnan(y) & x>1))));
                R2_VPD(k)=R(1,2)^2;
             
            else
                RMSE(k)=NaN;
                par_gs_VPD(k,:) = [NaN NaN];
                R2_VPD(k)=NaN;
            end    
           
            VPD_BOUND{k}=VPD_UP;
            GS_BOUND{k}=GS_UP;
            clear GS_UP VPD_UP
            
end
%% mulitple constraint function (radiation)
% SWIN: incoming shortwave radiation
for k=1:lng;
    k
    % load here variables for site
    % SWIN, VPD, DOY, HH, LE, H, GS, ga
    AV=LE+H; % available energy flux
    GS(GS<0)=NaN; % remove negative surface conductance
    
    % create indices for selection
    rem = find((DOY<DOY_start | DOY>=DOY_end) ...
            | (HH<time_start | HH>=time_end) | AV<25);
    
    VPD(rem)=NaN;
    GS(rem)=NaN;
    SWIN(rem)=NaN;
    
    % SWin bins for parameter fitting
    PRC_SWIN=50:100:1150; % 100 Wm-2 bins
    i=find(~isnan(GS));
    SWIN_UP = [];
    GS_UP = [];
    for l =1:length(PRC_SWIN)+1;
            if l == 1;
                ind=find(SWIN<=PRC_SWIN(l));     
            elseif l<=length(PRC_SWIN);
                ind=find(SWIN>PRC_SWIN(l-1) & SWIN<=PRC_SWIN(l));  
            else
                ind=find(SWIN>PRC_SWIN(l-1));
            end
            
            if length(ind)>3 & ~isempty(i)
            	
                VPD_sel = VPD(ind);
                GS_sel = GS(ind);
                SWIN_sel = SWIN(ind);
                
                % select outliers using quartile method (GS)
                IQR = prctile(GS_sel,75)-prctile(GS_sel,25);
                outl = find(GS_sel>prctile(GS_sel,75)+1.5*IQR | GS_sel<prctile(GS_sel,5)-1.25*IQR);
                GS_sel(outl)=NaN;

                SWIN_BINS(k,l)=nanmean(SWIN_sel);
                GS_BINS(k,l)=nanmean(GS_sel);
                GS_BINS_SD(k,l)=nanstd(GS_sel);
                
                % upper boundary
                x1=SWIN_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
                y1=GS_sel(GS_sel>=GS_BINS(k,l)+1.*GS_BINS_SD(k,l));
              
                clear GS_sel VPD_sel RN_sel
                SWIN_BINS_UP(k,l)=nanmedian(x1);
                GS_BINS_SW_UP(k,l)=nanmedian(y1);
                SWIN_UP=vertcat(SWIN_UP,x1);
                GS_UP=vertcat(GS_UP,y1);
                clear x1 y1 w1
            else
                SWIN_BINS_UP(k,l)=NaN;
                GS_BINS_SW_UP(k,l)=NaN;
            end
    end 

    % fit model for light-gs relationship
            if sum(~isnan(GS_UP))>10;
                x1 = SWIN_UP;
                y = (GS_UP./GS_max(k));
                i1 = find(~isnan(x1) & ~isnan(y));
                
                gsFxnSWIN=@(params,SWIN) ((params(1).*params(2).*SWIN)./(params(1)+params(2).*SWIN)); % where params is FCmax, alpha, RecoBase, Q10

                par0_gs=[15 0.05];
                par_gs_SWIN(k,:) = nlinfit(x1(i1),y(i1),gsFxnSWIN,par0_gs,opts);

                y=GS_BINS_SW_UP(k,:)./GS_max(k);
                x=SWIN_BINS_UP(k,:);
                RMSE_SWIN(k)=sqrt(mean(y(~isnan(x) & ~isnan(y) & x>1)-(gsFxnSWIN(par_gs_SWIN(k,:),x(~isnan(x) &~isnan(y) & x>1)))).^2);
                R=corrcoef(y(~isnan(x) & ~isnan(y) & x>1),(gsFxnSWIN(par_gs_SWIN(k,:),x(~isnan(x) &~isnan(y) & x>1))));
                R2_SWIN(k)=R(1,2)^2;
            else
                RMSE_SWIN(k)=NaN;
                par_gs_SWIN(k,:) = [NaN NaN];
                R2_SWIN(k)=NaN;
            end    
            
            SWIN_BOUND{k}=SWIN_UP;
            GS_BOUND{k}=GS_UP;
            clear GS_UP SWIN_UP

end
close all

%% test gs model
for k=1:lng;
    % load DOY, HH, VPD, GS, SWIN
    AV = H+LE;
    rem = find((DOY<DOY_start | DOY>=DOY_end)...
            | (HH<time_start | HH>=time_end) | AV<25);
    
    GS(GS<0)=NaN;
    
    VPD(rem)=NaN;
    GS(rem)=NaN;
    SWIN(rem)=NaN;
    
    if ~isnan(par_gs_VPD(k,1)) & ~isnan(par_gs_SWIN(k,1))
        x=GS_max(k).*(gsFxn(par_gs_VPD(k,:),VPD)).*(gsFxnSWIN(par_gs_SWIN(k,:),SWIN));
        x=real(x);
        % remove data > 99th percentile
        x(x>prctile(x,99))=NaN;
        y=GS;
        y=real(y);
        % remove data > 99th percentile
        y(y>prctile(y,99))=NaN;
        % calculate root-mean-squared-error
        RMSE_FIT(k)=sqrt(nanmean((x-y).^2));
        ind=find(~isnan(x) & ~isnan(y));
        % calculate R2
        [R p] = corrcoef(x(ind),y(ind));
        R2_FIT(k)=R(1,2)^2;
    else
        R2_FIT(k)=NaN;
        RMSE_FIT(k)=NaN;
    end
    % plot EC-derived gs vs. modelled gs
    figure,
    scatter(GS_max(k).*(gsFxn(par_gs_VPD(k,:),VPD)).*(gsFxnSWIN(par_gs_SWIN(k,:),SWIN)),GS,'.')
    hold on
    plot(0:0.0005:0.02,0:0.0005:0.02)
    xlim([0 0.05]);
    ylim([0 0.05]);
    title(num2str(k))
    pause(1)
    close all
end
R2_FIT(R2_FIT==0)=NaN;
