%Run after Cloud_Model_ECS
CloudEBM_model_ECS_Steady_State


%Coefficients_linear = polyfit(Annual_T_obs, epsilon_ClearSky_obs,1)

dR_dT_cloudless_ConstHoriz = 4*sigma*Coefficients_linear(2)*Annual_T3_obs + 5*sigma*Coefficients_linear(1)*Annual_T4_obs;

dT_dR_cloudless_ConstHoriz = 1./dR_dT_cloudless_ConstHoriz;
dT_dR_cloudless_ConstHoriz_equatorial = 0.5*(dT_dR_cloudless_ConstHoriz(18)+dT_dR_cloudless_ConstHoriz(19));




dR_dT_cloudless_wSurf_albedo = 4*sigma*Coefficients_linear(2)*Annual_T3_obs + 5*sigma*Coefficients_linear(1)*Annual_T4_obs + S_in_start.*function_dsurface_alphadT(Annual_T_obs, 36);


figure(1)
hold on;
plot(phi_deg, (1./(4*sigma*Annual_T3_obs))/(0.5*(1/(4*sigma*Annual_T3_obs(18) )+1/((4*sigma*Annual_T3_obs(19))))))
plot(phi_deg, dT_dR_cloudless_ConstHoriz/dT_dR_cloudless_ConstHoriz_equatorial)
X=[LGM_phi_deg', fliplr(LGM_phi_deg')];
Y=[Normalised_LGM_amp_low, fliplr(Normalised_LGM_amp_high)];
h2 = fill(X,Y, 'k', 'facealpha', 0.2);

plot (LGM_phi_deg, Normalised_LGM_amp, 'k', 'linewidth',2.0)
plot(phi_HadCRUT5, (Mean_2000_2020 - Mean_1900_1920)/(0.5*(Mean_2000_2020(18) + Mean_2000_2020(19) - Mean_1900_1920(18) -Mean_1900_1920(19) ) ) );




function dsurface_alphadT = function_dsurface_alphadT(T, size)

    dsurface_alphadT = 3.0*4.77471824e-6*T.*T - 2.0*3.73526338e-3*T + 9.60436190e-1 ;
    %Need to cycle through index
    for i=1:size
        if(dsurface_alphadT(i) > 0.0)
            dsurface_alphadT(i) = 0.0;
        end
    end

end
