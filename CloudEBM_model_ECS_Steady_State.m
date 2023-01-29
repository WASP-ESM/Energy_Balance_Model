%Energy Balance model with latitude-varying cloud albedo, temperature sensitive surface albedo, and temperature-sensitive emissivity (WVLR feedback).

%Define the turning points for the temperature-surface albedo link
T_cold=243.15;
T_warm=273.15-5.0;
alpha_mean_cold = 0.45;
alpha_mean_warm = 0.17;

%%%%%Loading real-world data
Load_Data_EBM;
%Load_HadCRUT5_LGM %Need to download HadCRUT5 data and
%Load_PETM %Need to download PETM data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_deg = linspace(-87.5, 87.5, (180/5));
phi_rad = (pi()/180)* phi_deg;

[~, size_array] = size(phi_rad)
phi_deg_between = linspace(0, 0, size_array -1);
for i=1:size_array-1
    phi_deg_between(i) = phi_deg(i)+0.5*180/size_array;
end


S_in_daily_calc=linspace(0, 0, size_array);
for i=1:12
    for j=1:size_array
        S_in_daily_calc(j) = S_in_daily_calc(j) + Insolation_monthly_average(j,i);
    end
end
S_in_from_daily_calc = S_in_daily_calc/12;

S_in = S_in_from_daily_calc;

S_in_start = S_in;


%Read in CRU absolute temperature record by monthly average
CRUTEM_Kelvin(:,:,1)=CRUTEM_Jan + 273.15;
CRUTEM_Kelvin(:,:,2)=CRUTEM_Feb + 273.15;
CRUTEM_Kelvin(:,:,3)=CRUTEM_Mar + 273.15;
CRUTEM_Kelvin(:,:,4)=CRUTEM_Apr + 273.15;
CRUTEM_Kelvin(:,:,5)=CRUTEM_May + 273.15;
CRUTEM_Kelvin(:,:,6)=CRUTEM_Jun + 273.15;
CRUTEM_Kelvin(:,:,7)=CRUTEM_Jul + 273.15;
CRUTEM_Kelvin(:,:,8)=CRUTEM_Aug + 273.15;
CRUTEM_Kelvin(:,:,9)=CRUTEM_Sep + 273.15;
CRUTEM_Kelvin(:,:,10)=CRUTEM_Oct + 273.15;
CRUTEM_Kelvin(:,:,11)=CRUTEM_Nov + 273.15;
CRUTEM_Kelvin(:,:,12)=CRUTEM_Dec + 273.15;

T4_CRUTEM = CRUTEM_Kelvin.^4;
T3_CRUTEM = CRUTEM_Kelvin.^3;

dt = 60*60*24; %time-step (1 'day', note there is no seasonality in the model)

R_Earth = 6371000.0 % Radius of Planet earth in m
Area_Earth = 4*pi()*R_Earth^2;



Length_phi = linspace(0,0,size_array-1);
for i=1:size_array-1
    Length_phi(i) = 2*pi()*R_Earth * cos(phi_rad(i)+0.5*pi()/size_array); %length of line of latitude mid way between the points
end

Area = 0.5*Area_Earth * ( (1-sin(phi_rad-0.5*pi()/size_array)) - (1-sin(phi_rad+0.5*pi()/size_array)) ); %Areas of latitudinal bands

sigma = 5.67e-8; %Stephan-Botlzmann constant in kg s^-3 K^-4

Annual_T4_obs = linspace(0, 0, size_array);
Annual_T_obs = linspace(0, 0, size_array);
Annual_T3_obs = linspace(0,0,size_array);
%Now take annual and latitudinal average of CRU T^4 etc

for i = 1:72
    for j = 1:12
        for k = 1:size_array
            Annual_T4_obs(k) = Annual_T4_obs(k) + (T4_CRUTEM(k,i,j));
            Annual_T_obs(k) = Annual_T_obs(k) + CRUTEM_Kelvin(k,i,j);
            Annual_T3_obs(k) = Annual_T3_obs(k) + T3_CRUTEM(k,i,j);
        end
    end
end
Annual_T4_obs = Annual_T4_obs/(12*72);
Annual_T_obs = Annual_T_obs/(12*72);
Annual_T3_obs = Annual_T3_obs/(12*72);

Annual_T_obs_interp = linspace(0,0,size_array-1);

for(i=1:size_array-1)
    Annual_T_obs_interp(i) = 0.5*(Annual_T_obs(i) + Annual_T_obs(i+1));
end


epsilon_AllSky_obs = L_out_obs./(sigma*Annual_T4_obs);
alpha_AllSky_obs = S_out_obs./S_in;

epsilon_ClearSky_obs = L_out_ClearSky_obs./(sigma*Annual_T4_obs);
alpha_ClearSky_obs = S_out_ClearSky_obs./S_in;

epsilon_CloudySky_obs = (1./Cloud_Amount).* (epsilon_AllSky_obs - (1-Cloud_Amount).*epsilon_ClearSky_obs );

c_epsilonCloudy = (1-epsilon_CloudySky_obs)./(1-epsilon_ClearSky_obs);

c_epsilonCloudy_obs = c_epsilonCloudy;

alpha_CloudySky_obs = (1./Cloud_insolation_fraction).*(alpha_AllSky_obs - (1-Cloud_insolation_fraction).*alpha_ClearSky_obs);


Coefficients_linear = polyfit(Annual_T_obs, epsilon_ClearSky_obs,1);


alpha_Cloud_obs = linspace(0.0,0.0,size_array);

Mean_c_epsilonCloudy = sum(c_epsilonCloudy.*Area)/sum(Area);

%Set heat capacities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 3.910e3; %Ocean mean heat capacity in J kg-1 K-1 (Williams et al., 2012)

Mass_ocean = 1027.0*sum(Area)*0.7*(2400-150); %Gregory (2000) depth of 2400m for lower ocean (ignores deep isolated basins) %1.35e21; % actual mass of ocean in kg

c_per_m2 = 0.7*150*c_p*1027.0; % Heat capacity of surface per metre squared (atmosphere plus ocean surface mixed layer) - 70% ocean coverage with 150m deep surface-ocean mixed layer. Ocean mixed layer dominated.

c = c_per_m2 * Area; %Heat capacity of each latitudinal band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Cloud_insolation_monthly=Insolation_monthly_average.*Cloud_Amount_AVHRR_5degree;




%%If using uniform cloud amount and cloud insolation fraction
Cloud_Amount = linspace(0, 0, size_array);
Cloud_insolation_fraction = linspace(0, 0, size_array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:12
    for j=1:size_array
        Cloud_Amount(j)=Cloud_Amount(j) + Cloud_Amount_AVHRR_5degree(j,i);
        Cloud_insolation_fraction(j) = Cloud_insolation_fraction(j) + Cloud_insolation_monthly(j,i);
    end
end

Cloud_Amount = Cloud_Amount/12;
Cloud_insolation_fraction = Cloud_insolation_fraction./(12*S_in);

Cloud_Amount_obs = Cloud_Amount;
Cloud_insolation_fraction_obs = Cloud_insolation_fraction;


%%Set diffusivity

kappa_poleward = linspace(90.0, 90.0, size_array-1); %horizontal poleward heat transport diffusivity in W per °C per m^2


alpha = linspace(0.3,0.3,size_array);

%declare initial surface temperature and northward heat trasport%%%%%%%%%%%%%
T_surface = linspace(273.15, 273.15, size_array);
Northward_Heat_Transport = linspace(0.0,0.0, size_array-1);

%Set initial clear sky emissivity, allsky emissivity, outgoing longwave and outgoing shortwave %%%%%%%%
epsilon_Clear_Sky_calc = 1.82785897 - 3.95009011e-03*T_surface;
epsilon_Cloudy_calc = 1-c_epsilonCloudy.*(1-epsilon_Clear_Sky_calc);
epsilon = Cloud_Amount.*epsilon_Cloudy_calc + (1-Cloud_Amount).*epsilon_Clear_Sky_calc; %1-1.3*(1-epsilon_Clear_Sky_calc);
L_out = sigma*epsilon.*(T_surface.^4);
S_out =alpha.*S_in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_epsilon = sum(epsilon_AllSky_obs.*Area)/sum(Area);

%%Set diffusivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa_poleward_obs = linspace(0,0,size_array-1);
for i=1:size_array-1
    kappa_poleward_obs(i) = - Northward_heat_flux_obs(5*i) /  ( (pi()*R_Earth/size_array)*(Annual_T_obs(i+1)-Annual_T_obs(i))*Length_phi(i) );
end

kappa_poleward = linspace(100.0, 100.0, size_array-1); %horizontal poleward heat transport diffusivity in W per °C per m^2



for i=1:size_array-1
    if(kappa_poleward_obs(i) > 0.0)
        kappa_poleward(i) = kappa_poleward_obs(i);
    end
    if(kappa_poleward_obs(i) < 0.0)
        kappa_poleward(i) = kappa_poleward(i-1);
    end
    %if(kappa_poleward(i) > 200.0)
    %    kappa_poleward(i) = 200.0;
    %end

end
kappa_standard = kappa_poleward;
kappa_excluding_tropics = [kappa_standard(1:16), kappa_standard(22:35)];
T_exluding_tropics = [Annual_T_obs(1:16), Annual_T_obs(22:35)];
Coeff_kappa_T_extratropics = polyfit(T_exluding_tropics,kappa_excluding_tropics, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%definition of radiative forcing
RF = linspace(0.0, 0.0, size_array);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 2.5e6; %J/kg latent heat water vapour
R_v = 461.0; %gas constant for water vapour in J per K per kg
T = linspace(200,320); %
c_p_dryair = 3910; %specific heat capacity of air J per K per kg
H_rel = 0.7;

qstar = 0.6*6.11*exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp)));

dqstar_dT = ((0.6*6.11*L)./(R_v.*Annual_T_obs_interp.*Annual_T_obs_interp)) .* exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp)));

r_latent_dry = (H_rel*L)/(1000*c_p_dryair) *dqstar_dT;

dr_ld_dT = (H_rel*L)/(1000*c_p_dryair) * (0.6*6.11*L)./(R_v*Annual_T_obs_interp.*Annual_T_obs_interp.*Annual_T_obs_interp).*exp((L/R_v)*((1/273) - (1./Annual_T_obs_interp))).*((L./(R_v*Annual_T_obs_interp)) - 2);

kappa_init = kappa_poleward_obs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



attempt = 0.5;
for i=1:size_array
    for j=1:100 %n-step iteration
        alpha_cloud_attempt = attempt + alpha_ClearSky_obs(i)*(1-2*attempt + attempt*attempt);
        attempt = attempt + (alpha_CloudySky_obs(i) - alpha_cloud_attempt)/2.0;
    end
    alpha_Cloud_obs(i) = attempt;
    attempt = 0.5;
end




%%%Contact with deep ocean%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_phi = linspace(0.0,0.0,size_array); %Wm-2K-1 from upwelling water
gamma = 1.6 ; %Wm-2K-1 global mean%%%%%%%
Area_upwelling = 0.0;
for(i=1:size_array)
    if(phi_deg(i) < -50.0 & phi_deg(i) > -70.0)
        Area_upwelling = Area_upwelling + Area(i);
    end
end

%62.5 per cent of sub-surface water next makes contact with surface mixed layer in Southern Ocean, 50 to 70 ° South (DeVries and Primeau, 2011)
gamma_upwellSO = 0.625*gamma*sum(Area)/Area_upwelling;
gamma_upwell_rest = (1-0.625)*gamma*sum(Area)/(sum(Area)-Area_upwelling);


DeltaT_deepocean = 0.0;

for i=1:size_array
    if (phi_deg(i) < -50.0 & phi_deg(i) > -70.0)
        gamma_phi(i) = gamma_upwellSO;
    else
        gamma_phi(i) = gamma_upwell_rest;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_Cloud_obs_planetary = (sum(Area.*alpha_Cloud_obs)/sum(Area));

%%Run model to equilibrium%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for model_type =1:7


for forcing = 1:2
    if(forcing == 1)
        T_surface = linspace(273.15, 273.15, size_array);
        Northward_Heat_Transport = linspace(0.0,0.0, size_array-1);

        %%Reset clouds
        Cloud_Amount = Cloud_Amount_obs;
        Cloud_insolation_fraction = Cloud_insolation_fraction_obs;
        
        %c_epsilonCloudy = c_epsilonCloudy_obs;
        for i=1:size_array
            c_epsilonCloudy(i) = sum(Area.*c_epsilonCloudy_obs)/sum(Area);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
   
    
    DeltaT_deepocean = 0.0;

    %%%%%%%%%%%%%%%
    tmax = 365*100; %initial spin up (number of 'days' of interation)

    if(model_type == 4 & forcing == 2)
        tmax = 365*20;
    end
    if(model_type == 5 & forcing == 2)
        tmax = 365*100;
    end
    if(model_type == 6 & forcing == 2)
        tmax = 365*500
    end
    if(model_type == 7 & forcing == 2)
        tmax = 365*1000
    end

    %%%%%%%%%%%%%%%%%

    for t=1:tmax %time loop, number of days

        

        if(forcing == 2 & model_type >= 1)
            for i=1:size_array
                c_epsilonCloudy(i) = sum(Area.*c_epsilonCloudy_obs)/sum(Area);
            end
        end
        
        if(forcing == 1)
            alpha_cloud = alpha_Cloud_obs_planetary*(1.0 + (1-alpha_Cloud_obs_planetary)*(0.5*(3*sin(phi_rad).*sin(phi_rad) - 1.0)));  %
        end
        

        if(forcing == 1)
            kappa_poleward = kappa_standard;
        end
        if(forcing == 2 & model_type>=2)
            
            for i=1:size_array-1
                if(kappa_standard(i) > 200) %do not apply over tropics
                    kappa_poleward(i) = kappa_standard(i);
                end
                if(kappa_standard(i) <= 200)
                    kappa_poleward(i) = kappa_standard(i) * (1 + function_ratio_Latent_Dry(T_surface(i), T_surface(i)+1.0, 0.70)) / (1 + function_ratio_Latent_Dry(T_surface1(i), T_surface1(i)+1.0, 0.70));
                    
                end
            end
        end

        %Transient deep ocean exchanges%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Heating_deep = 0.0;
        Heating_surface = linspace(0.0,0.0,size_array);
        if(forcing == 2 & model_type >=4 )
            for i=1:size_array
                Heating_surface(i) = - gamma_phi(i)*(T_surface(i) - T_spinup(i)-DeltaT_deepocean);
                Heating_deep = Heating_deep + dt*(Area(i)*gamma_phi(i)*(T_surface(i) - T_spinup(i)-DeltaT_deepocean));
            end
        end
        DeltaT_deepocean = DeltaT_deepocean + Heating_deep/(c_p*Mass_ocean);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %calculate horizontal heat transport%%%%%%%%%%%%%%%%%%%%%%
        for i=1:size_array-1
            Northward_Heat_Transport(i) = -(kappa_poleward(i)*pi()*R_Earth/size_array)*(T_surface(i+1)-T_surface(i))*Length_phi(i);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        %Update temperatures %%%%%%%%%
        T_surface(1) = T_surface(1) + dt*(1/c(1))*( Area(1)*(S_in(1) + RF(1) - L_out(1) - S_out(1) + Heating_surface(1)) - Northward_Heat_Transport(1));
        for i=2:size_array-1
            T_surface(i) = T_surface(i) + dt*(1/c(i))*(Area(i)*(S_in(i) + RF(i) - L_out(i) - S_out(i) + Heating_surface(i)) + Northward_Heat_Transport(i-1) - Northward_Heat_Transport(i));
        end
        T_surface(size_array) = T_surface(size_array) + dt*(1/c(size_array))*(Area(size_array)*(S_in(size_array) + RF(size_array) - L_out(size_array) - S_out(size_array) + Heating_surface(size_array)) + Northward_Heat_Transport(size_array-1));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Update albedo and emissivity %%%%%%%%%%%%%%%
    
        for i=1:size_array

            
            
            alpha_surface(i) =function_surface_alpha_phi(phi_rad(i), function_mean_alpha_T(T_surface(i), T_cold,alpha_mean_cold, T_warm,alpha_mean_warm));

    
            alpha(i) = alpha_surface(i)*(1.0-Cloud_insolation_fraction(i)) + Cloud_insolation_fraction(i)*(alpha_cloud(i) + alpha_surface(i)*(1.0 - alpha_cloud(i))*(1.0 - alpha_cloud(i)) );
        

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%Calculate emissivity in relation to surface temperature%%%%%%%%%%
        epsilon_Clear_Sky_calc = function_epsilon(T_surface, Coefficients_linear(1),Coefficients_linear(2)) ;


        for i=1:size_array
            if(epsilon_Clear_Sky_calc(i) > 0.96)
                epsilon_Clear_Sky_calc(i) = 0.96;
            end
            if(epsilon_Clear_Sky_calc(i) < 0.4)
            epsilon_Clear_Sky_calc(i) = 0.4;
            end
            epsilon_Cloudy_calc(i) = 1.0 - c_epsilonCloudy(i)*(1-epsilon_Clear_Sky_calc(i));
            %Cloud impact
            epsilon(i) = (1-Cloud_Amount(i))*epsilon_Clear_Sky_calc(i) + Cloud_Amount(i)*epsilon_Cloudy_calc(i);


        end
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Update vertical energy balance%%
    
        S_out =alpha.*S_in_start;
        L_out = sigma*epsilon.*(T_surface.^4);

    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Radiative forcing
        a_CO2_above285 = 5.35;%
        S_in = S_in_from_daily_calc;
        if (forcing == 1)
            CO2_logchange = log(1);
        end
        if (forcing == 2)
            if(model_type == 1)
                CO2_logchange = log(2.0);
            end
            if(model_type == 2)
                CO2_logchange = log(2.0);
            end
            if(model_type == 3)
                CO2_logchange = log(2.0);
            end
            if(model_type == 4)
                CO2_logchange = log(2.0);
            end
            if(model_type >= 5)
                CO2_logchange = log(2.0);
            end
        end
        for i=1:size_array
            RF(i) =  a_CO2_above285*CO2_logchange;
        end
        
        %Calculate northward heat transport per metre
        f_NHT = Northward_Heat_Transport./Length_phi;
        
        %calc df_NHT/dy
        df_NHT_dy = linspace(0.0,0.0, size_array);
        for i=2:size_array-1
            df_NHT_dy(i) = (Northward_Heat_Transport(i-1)-Northward_Heat_Transport(i))/(2*pi()*R_Earth*cos(phi_rad(i))*R_Earth*(phi_rad(i)-phi_rad(i-1))); %(f_NHT(i-1)-f_NHT(i))/(R_Earth*(phi_rad(i)-phi_rad(i-1)));
        end
        df_NHT_dy(1) = (-Northward_Heat_Transport(1))/(2*pi()*R_Earth*cos(phi_rad(1))*R_Earth*(phi_rad(2)-phi_rad(1)));
        df_NHT_dy(size_array) = (Northward_Heat_Transport(size_array-1))/(2*pi()*R_Earth*cos(phi_rad(size_array))*R_Earth*(phi_rad(size_array)-phi_rad(size_array-1))); f_NHT(size_array-1)/(R_Earth*(phi_rad(2)-phi_rad(1)));
        
        
        
    end %End of loop for timesteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(forcing == 1)
        T_spinup = T_surface;
    end
    
    if(model_type == 1 & forcing == 1)
        T_surface1 = T_surface;
        L_out1 = L_out;
        S_out1 = S_out;
        Northward_Heat_Transport1 = Northward_Heat_Transport;
        alpha1 = alpha;
        epsilon1 = epsilon;
        f_NHT1 = f_NHT;
        df_NHT_dy1 = df_NHT_dy;
        kappa1 = kappa_poleward;
        Heating_surface1 = Heating_surface;
    end

    if(model_type == 2 & forcing == 1)
        T_surface2 = T_surface;
        L_out2 = L_out;
        S_out2 = S_out;
        Northward_Heat_Transport2 = Northward_Heat_Transport;
        alpha2 = alpha;
        epsilon2 = epsilon;
        f_NHT2 = f_NHT;
        df_NHT_dy2 = df_NHT_dy;
        kappa2 = kappa_poleward;
        Heating_surface2 = Heating_surface;
    end

    if(model_type == 3 & forcing == 1)
        T_surface3 = T_surface;
        L_out3 = L_out;
        S_out3 = S_out;
        Northward_Heat_Transport3 = Northward_Heat_Transport;
        alpha3 = alpha;
        epsilon3 = epsilon;
        f_NHT3 = f_NHT;
        df_NHT_dy3 = df_NHT_dy;
        kappa3 = kappa_poleward;
        Heating_surface3 = Heating_surface;
    end

    if(model_type == 4 & forcing == 1)
        T_surface4 = T_surface;
        L_out4 = L_out;
        S_out4 = S_out;
        Northward_Heat_Transport4 = Northward_Heat_Transport;
        alpha4 = alpha;
        epsilon4 = epsilon;
        f_NHT4 = f_NHT;
        df_NHT_dy4 = df_NHT_dy;
        kappa4 = kappa_poleward;
        Heating_surface4 = Heating_surface;
    end

    if(model_type == 5 & forcing == 1)
        T_surface5 = T_surface;
        L_out5 = L_out;
        S_out5 = S_out;
        Northward_Heat_Transport5 = Northward_Heat_Transport;
        alpha5 = alpha;
        epsilon5 = epsilon;
        f_NHT5 = f_NHT;
        df_NHT_dy5 = df_NHT_dy;
        kappa5 = kappa_poleward;
        Heating_surface5 = Heating_surface;
    end

    if(model_type == 6 & forcing == 1)
        T_surface6 = T_surface;
        L_out6 = L_out;
        S_out6 = S_out;
        Northward_Heat_Transport6 = Northward_Heat_Transport;
        alpha6 = alpha;
        epsilon6 = epsilon;
        f_NHT6 = f_NHT;
        df_NHT_dy6 = df_NHT_dy;
        kappa6 = kappa_poleward;
        Heating_surface6 = Heating_surface;
    end

    if(model_type == 7 & forcing == 1)
        T_surface7 = T_surface;
        L_out7 = L_out;
        S_out7 = S_out;
        Northward_Heat_Transport7 = Northward_Heat_Transport;
        alpha7 = alpha;
        epsilon7 = epsilon;
        f_NHT7 = f_NHT;
        df_NHT_dy7 = df_NHT_dy;
        kappa7 = kappa_poleward;
        Heating_surface7 = Heating_surface;
    end

if(model_type == 1 & forcing == 2)
    T_surface8 = T_surface;
    L_out8 = L_out;
    S_out8 = S_out;
    Northward_Heat_Transport8 = Northward_Heat_Transport;
    alpha8 = alpha;
    epsilon8 = epsilon;
    f_NHT8 = f_NHT;
    df_NHT_dy8 = df_NHT_dy;
    kappa8 = kappa_poleward;
    Heating_surface8 = Heating_surface;
end

if(model_type == 2 & forcing == 2)
    T_surface9 = T_surface;
    L_out9 = L_out;
    S_out9 = S_out;
    Northward_Heat_Transport9 = Northward_Heat_Transport;
    alpha9 = alpha;
    epsilon9 = epsilon;
    f_NHT9 = f_NHT;
    df_NHT_dy9 = df_NHT_dy;
    kappa9 = kappa_poleward;
    Heating_surface9 = Heating_surface;
end

if(model_type == 3 & forcing == 2)
    T_surface10 = T_surface;
    L_out10 = L_out;
    S_out10 = S_out;
    Northward_Heat_Transport10 = Northward_Heat_Transport;
    alpha10 = alpha;
    epsilon10 = epsilon;
    f_NHT10 = f_NHT;
    df_NHT_dy10 = df_NHT_dy;
    kappa10 = kappa_poleward;
    Heating_surface10 = Heating_surface;
end

if(model_type == 4 & forcing == 2)
    T_surface11 = T_surface;
    L_out11 = L_out;
    S_out11 = S_out;
    Northward_Heat_Transport11 = Northward_Heat_Transport;
    alpha11 = alpha;
    epsilon11 = epsilon;
    f_NHT11 = f_NHT;
    df_NHT_dy11 = df_NHT_dy;
    kappa11 = kappa_poleward;
    Heating_surface11 = Heating_surface;
end

if(model_type == 5 & forcing == 2)
    T_surface12 = T_surface;
    L_out12 = L_out;
    S_out12 = S_out;
    Northward_Heat_Transport12 = Northward_Heat_Transport;
    alpha12 = alpha;
    epsilon12 = epsilon;
    f_NHT12 = f_NHT;
    df_NHT_dy12 = df_NHT_dy;
    kappa12 = kappa_poleward;
    Heating_surface12 = Heating_surface;
end

if(model_type == 6 & forcing == 2)
    T_surface13 = T_surface;
    L_out13 = L_out;
    S_out13 = S_out;
    Northward_Heat_Transport13 = Northward_Heat_Transport;
    alpha13 = alpha;
    epsilon13 = epsilon;
    f_NHT13 = f_NHT;
    df_NHT_dy13 = df_NHT_dy;
    kappa13 = kappa_poleward;
    Heating_surface13 = Heating_surface;
end

if(model_type == 7 & forcing == 2)
    T_surface14 = T_surface;
    L_out14 = L_out;
    S_out14 = S_out;
    Northward_Heat_Transport14 = Northward_Heat_Transport;
    alpha14 = alpha;
    epsilon14 = epsilon;
    f_NHT14 = f_NHT;
    df_NHT_dy14 = df_NHT_dy;
    kappa14 = kappa_poleward;
    Heating_surface14 = Heating_surface;
end




end

end


figure(6)
subplot(3,2,1)
plot(phi_deg, S_in-L_out1-S_out1, phi_deg, S_in-L_out_obs-S_out_obs )
subplot(3,2,2)
plot(phi_deg, L_out1, phi_deg, S_out1, phi_deg, L_out_obs, phi_deg, S_out_obs)
subplot(3,2,3)
plot(phi_deg, T_surface1, phi_deg, Annual_T_obs)
subplot(3,2,4)
plot(phi_deg_between, Northward_Heat_Transport1, phi_1deg, Northward_heat_flux_obs)
subplot(3,2,5)
plot(phi_deg, alpha1, phi_deg, alpha_AllSky_obs)
subplot(3,2,6)
plot(phi_deg, epsilon1, phi_deg, epsilon_AllSky_obs)

    

T_range = linspace(min(Annual_T_obs), max(Annual_T_obs), 500);

figure(3)
plot(Annual_T_obs, epsilon_ClearSky_obs, T_range, function_epsilon(T_range, Coefficients_linear(1),Coefficients_linear(2)))


figure(2)
subplot(4,2,1)
plot(phi_deg, S_in_start, phi_deg, S_out_obs, phi_deg, S_out_ClearSky_obs)
subplot(4,2,2)
plot(phi_deg, sigma*Annual_T4_obs , phi_deg, L_out_obs, phi_deg, L_out_ClearSky_obs)
subplot(4,2,3)
plot(phi_deg, epsilon_ClearSky_obs, phi_deg, epsilon_AllSky_obs, phi_deg, epsilon_CloudySky_obs)
subplot(4,2,4)
plot(phi_deg, alpha_ClearSky_obs, phi_deg, alpha_AllSky_obs, phi_deg, alpha_CloudySky_obs)
subplot(4,2,7)
plot(phi_deg, c_epsilonCloudy_obs, phi_deg, c_epsilonCloudy)
subplot(4,2,5)
plot(phi_deg, Cloud_Amount_obs )
subplot(4,2,6)
plot(phi_deg, Cloud_Amount_obs, phi_deg, Cloud_insolation_fraction_obs)
subplot(4,2,8)
plot(phi_deg, alpha_Cloud_obs,phi_deg, alpha_cloud)




%For figure 4

alpha_phi_mean0=linspace(0.0,0.0,180);
alpha_phi_mean01=linspace(0.0,0.0,180);
alpha_phi_mean02=linspace(0.0,0.0,180);
alpha_phi_mean03=linspace(0.0,0.0,180);
alpha_phi_mean04=linspace(0.0,0.0,180);
alpha_phi_mean05=linspace(0.0,0.0,180);
alpha_phi_mean06=linspace(0.0,0.0,180);
alpha_phi_mean07=linspace(0.0,0.0,180);
alpha_phi_mean08=linspace(0.0,0.0,180);
alpha_phi_mean09=linspace(0.0,0.0,180);
alpha_phi_mean1=linspace(0.0,0.0,180);

T_fig12 = linspace(220, 320, 100);
alpha_mean_fig12 = linspace(0.0,0.0,100);

for(i=1:180)
    alpha_phi_mean0(i) = function_surface_alpha_phi(deg2rad(phi_1deg   (i)),0.0);
    alpha_phi_mean01(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.1);
    alpha_phi_mean02(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.2);
    alpha_phi_mean03(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.3);
    alpha_phi_mean04(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.4);
    alpha_phi_mean05(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.5);
    alpha_phi_mean06(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.6);
    alpha_phi_mean07(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.7);
    alpha_phi_mean08(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.8);
    alpha_phi_mean09(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),0.9);
    alpha_phi_mean1(i) = function_surface_alpha_phi(deg2rad(phi_1deg(i)),1.0);
end

for(i=1:100)
    alpha_mean_fig12(i) = function_mean_alpha_T(T_fig12(i), T_cold,alpha_mean_cold, T_warm,alpha_mean_warm);
end

figure(4)
subplot(2,1,1)
plot(T_fig12, alpha_mean_fig12, T_cold, alpha_mean_cold, T_warm,alpha_mean_warm);
subplot(2,1,2)
plot(phi_1deg, alpha_phi_mean1, phi_1deg, alpha_phi_mean09, phi_1deg, alpha_phi_mean08, phi_1deg, alpha_phi_mean07, phi_1deg, alpha_phi_mean06, phi_1deg, alpha_phi_mean05, phi_1deg, alpha_phi_mean04, phi_1deg, alpha_phi_mean03, phi_1deg, alpha_phi_mean02, phi_1deg, alpha_phi_mean01, phi_1deg, alpha_phi_mean0, phi_deg, alpha_Cloud_obs);


figure(5)
hold on;

subplot(3,2,1)
plot(phi_deg_between, kappa8, phi_deg_between, kappa10, phi_deg_between, kappa11, phi_deg_between, kappa_poleward_obs)

subplot(3,2,2)
hold on;
%plot(phi_HadCRUT5, (Mean_2000_2020 - Mean_1900_1920)/(0.5*(Mean_2000_2020(18) + Mean_2000_2020(19) - Mean_1900_1920(18) -Mean_1900_1920(19) ) ) );
plot( phi_deg, (T_surface8-T_surface1)/(0.5*(T_surface8(18)+T_surface8(19)-T_surface1(18)-T_surface1(19))), phi_deg, (T_surface10-T_surface3)/(0.5*(T_surface10(18)+T_surface10(19)-T_surface3(18)-T_surface3(19))), phi_deg, (T_surface11-T_surface4)/(0.5*(T_surface11(18)+T_surface11(19)-T_surface4(18)-T_surface4(19))));

subplot(3,2,3)
plot( phi_deg, T_surface8-T_surface1, phi_deg, T_surface10-T_surface3, phi_deg, T_surface11-T_surface4);

subplot(3,2,4)
plot(phi_deg, (S_out8+L_out8-S_out1-L_out1)./(T_surface8-T_surface1), phi_deg, (S_out10+L_out10-S_out3-L_out3)./(T_surface10-T_surface3), phi_deg, (S_out11+L_out11-S_out4-L_out4)./(T_surface11-T_surface4));

subplot(3,2,6)
plot(phi_deg, L_out8+S_out8 - RF, phi_deg, L_out10+S_out10 - RF, phi_deg, L_out11+S_out11 - RF, phi_deg, L_out1+S_out1)

subplot(3,2,5)
plot(phi_deg_between, f_NHT8, phi_deg_between, f_NHT10, phi_deg_between, f_NHT11, phi_deg_between, f_NHT1)


figure(9)
subplot(3,1,1)
plot(phi_deg, gamma_phi)
subplot(3,1,2)
plot(phi_deg, T_surface11-T_surface4, phi_deg, T_surface12-T_surface5, phi_deg, T_surface13-T_surface6, phi_deg, T_surface14-T_surface5)
subplot(3,1,3)
plot(phi_deg, Heating_surface11, phi_deg, Heating_surface12, phi_deg, Heating_surface13, phi_deg, Heating_surface14)


%functions

function f_ratio_Latent_Dry = function_ratio_Latent_Dry(T1, T2, RH)
    L_v = 2.5e6; %J/kg latent heat water vapour
    R_v = 461; %Gas constant for water vapour
    c_p_dryair = 1000.5; %specific heat capacity dry air at constant pressure in J/(kgK)
    qstar1 = (30/50)*6.11*exp(L_v/R_v*((1/273.0)- (1/T1)));
    qstar2 = (30/50)*6.11*exp(L_v/R_v*((1/273.0)- (1/T2)));
    Latent = RH*L_v*(qstar2-qstar1)/1000;
    Dry = c_p_dryair*(T2 - T1);
    f_ratio_Latent_Dry = Latent/Dry;
end

function f_epsilonClearSky = function_epsilon(T, A, B)
    f_epsilonClearSky = A*T + B;
end



function mean_alpha = function_mean_alpha_T(Temp, a,b,c,d)
    %define turning points (a,b) and (c,d)
    
    
    k = -6.0*(b-d)/((a-c)*(a-c)*(a-c));
    h = ((b+d)/2.0) - (( (b-d)*(a+c)*(a*a + c*c - 4.0*a*c) ) / (2.0*(a-c)*(a-c)*(a-c)) );
    mean_alpha = d;
    
    if(Temp > c)
        mean_alpha = d;
    end
    if(Temp > a & Temp <= c)
        mean_alpha = k*((Temp*Temp*Temp)/3.0 - (a+c)*Temp*Temp/2 + a*c*Temp) + h;
    end
    if(Temp <= a)
        mean_alpha = b;
    end

end



function surface_alpha = function_surface_alpha_phi(phi,  mean_alpha)

    surface_alpha = mean_alpha*(1.0 + (1.0-mean_alpha)*(0.5*(3.0*sin (phi) * sin (phi) - 1.0)));

end





