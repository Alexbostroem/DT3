clc
clear

%% D1
%Data given
Alti = 35000 * 0.3048; 
M = 0.82;
dTisa = 10;
BPR = 11.8;
FPR = 1.55;
OPR = 47;
HPC = 4.5;
T_inlet = 1680;
inlet_eff = 0.995;
fan_poly_eff = 0.92;
cc_eff = 0.999;
cc_pressure_loss = 0.035;
HPT_poly_eff = 0.905;
IPT_poly_eff = 0.91;
LPT_poly_eff = 0.915;
Cold_jet_eff = 0.98;
Hot_jet_eff = 0.99;
Shaft_mec_eff = 0.995;
IPC_poly_eff = 0.91;
HPC_poly_eff = 0.915;
Net_thrust = 67.3 * 10^3;


%Data from charts
TSL = 288.15;
PSL = 101325;
g = 9.81;
R = 287.05;
gamma_air = 1.4;
gamma_gas = 1.333;
cp_air = 1005;
cp_gas = 1148;



psl = 101325;
tsl = 288.15;
P_a = psl*(1-0.0065*Alti/tsl)^5.2561;
T_a = tsl-6.5*Alti/1000 +dTisa;


%Calculating stagnating at intake T01 AND P01 (station 1)
T_01 = T_a * (1 + ((gamma_air-1)/2)*M^2); %(1.30B FS)
P_01 = P_a * (1 + inlet_eff*((gamma_air-1)/2)*M^2)^(gamma_air/(gamma_air-1));

%Calculating station 2 (After Fan)
T_02 = T_01*(FPR^((gamma_air-1)/(fan_poly_eff*gamma_air))); %(1.32 FS)
P_02 = P_01*FPR;


%Calculating station 3 (After IPC)
IPR = OPR/(FPR*HPC);
T_03 = (IPR^((gamma_air-1)/(IPC_poly_eff*gamma_air)))*T_02; %(1.32 FS)
P_03 = P_02*IPR;


%Calculating station 4 (After HPC)
T_04 = (HPC^((gamma_air-1)/(HPC_poly_eff*gamma_air)))*T_03; %(1.32 FS)
P_04 = P_03*HPC;


%Calculation FAR (station CC)
FAR1 =0.10118 + 2.00376E-05*(700-T_04);
FAR2 =3.7078E-03 - 5.2368E-06*(700-T_04) - 5.2632E-06*T_inlet;
FAR3 =8.889E-08*(T_inlet-950);
FAR_real = (FAR1-sqrt(FAR1^2 + FAR2)-FAR3)/cc_eff;

m_a = 670.5154;

[m_cold,m_h,m_fuel,m_cooling_s,m_cooling_r,m_h_corr] = massflow_new(m_a,FAR_real,BPR,0.2);

%Calcualtion station 5 (After compustion)
T_05 = T_inlet;
P_05 = P_04 * (1 - cc_pressure_loss);


%Calculating station 6 (After HPT)

    %Calculating after stator (b)
    m_6a = m_h_corr + m_fuel;
    m_6b = m_6a + m_cooling_s;
    T_06b = ((m_6a*cp_gas*T_05) + (m_cooling_s*cp_air*T_04))/(m_6b*cp_gas);
    P_06b = P_05;

    %Calculating after rotor (c)
    T_06c = T_06b - ((m_h/m_6b) * (cp_air/cp_gas) * ((T_04-T_03)/Shaft_mec_eff));
    P_06c = P_06b/((T_06b/T_06c)^(gamma_gas/((gamma_gas-1)*HPT_poly_eff)));

    %Calculating mixing flow
    m_6 = m_6b + m_cooling_r;

 T_06 = ((m_6b*cp_gas*T_06c) + (m_cooling_r*cp_air*T_04))/(m_6*cp_gas); %Hade T0_5 för temp på rotor luft, borde va T_04?
 P_06 = P_06c;

%Calculating station 7 (After IPT)
T_07 = T_06 - (m_h/m_6 * (cp_air/cp_gas) * ((T_03-T_02)/Shaft_mec_eff));
P_07 = P_06/((T_06/T_07)^(gamma_gas/((gamma_gas-1)*IPT_poly_eff)));

%Calculating station 8 (After LPT)
T_08 = T_07 - (m_a/m_6 * (cp_air/cp_gas) * ((T_02-T_01)/Shaft_mec_eff));
P_08 = P_07/(T_07/T_08)^(gamma_gas/((gamma_gas-1)*LPT_poly_eff));

%Calculating stattion 9 (At nozzel outlet)
    %Calculating P_critical_hot(subtask)
    P_9_ratio_hot = P_08/P_a;
    P_9_critical_hot = 1/((1-(1/Hot_jet_eff)*((gamma_gas-1)/(gamma_gas+1)))^(gamma_gas/(gamma_gas-1)));

    %Calculating P_critical_cold(subtask)
    P_10_ratio_cold = P_02/P_a;
    P_10_critical_cold = 1/((1-((1/Cold_jet_eff)*((gamma_air-1)/(gamma_air+1))))^(gamma_air/(gamma_air-1)));
    
    if P_9_ratio_hot > P_9_critical_hot 
        %From (2.7) (Shocked)
        T_9 = 2*T_08/(gamma_gas+1);
        P_9 = P_08/P_9_critical_hot;
        C_9 = sqrt(gamma_gas*R*T_9);
        rho_9 = P_9 / (R*T_9);
        A_9 = (m_h + m_fuel) / (rho_9 * C_9);

        %Calculating thrust for hot flow
        FGH = (m_6) * C_9 + (A_9*(P_9 - P_a));
    else
        P_9 = P_a; %ÄNDRADE HÄR
        T_9 = T_08 - (Hot_jet_eff*T_08*(1 - (1/(P_08/P_9))^((gamma_gas-1)/gamma_gas)));
        C_9 = sqrt(2*cp_gas*(T_08-T_9));
      
        FGH = (m_6) * C_9;
    end

    if P_10_ratio_cold > P_10_critical_cold 
        %From (2.7) (Shocked)
        T_10 = 2*T_02/(gamma_air+1);
        P_10 = P_02/P_10_critical_cold;
        C_10 = sqrt(gamma_air*R*T_10);
        rho_10 = P_10 / (R*T_10);
        A_10 = m_cold / (rho_10 * C_10);
        FGC = m_cold * C_10 + (A_10*(P_10 - P_a));

      else
        P_10 = P_a;
        T_10 = T_02 - (Cold_jet_eff*T_02*(1 - (1/(P_02/P_10)^((gamma_air-1)/gamma_air))));
        C_10 = sqrt(2*cp_air*(T_02-T_10));
        FGC = m_cold * C_10;

    end

     C_a = M*sqrt(R*T_a*gamma_air); 
     FD = C_a*m_a;
     FNET = FGH + FGC - FD;



%% D2
    %% FAN
    % FAN Area
    M_ax = 0.603;
    fanArea = calculateFanArea(M_ax ,gamma_air,T_01,P_01,m_a);
    
    
    %FAN Rhub and Rtip 1
    EIS = 2020;
    v_fan = 44.29/(98.94+exp(0.01850*EIS-33.31));
    Rfan_tip1 = sqrt(fanArea/(pi*(1-v_fan^2)));
    Rfan_hub1 = v_fan * Rfan_tip1;

    %FAN Rhub and Rtip 2
    AR_fan = 2.4;
    alfa = 15;
    l_ax_fan = (1.98*Rfan_tip1 - (2*Rfan_hub1))/((2 * AR_fan) + tand(alfa));
    Rfan_hub2 = Rfan_hub1 + (l_ax_fan*tand(alfa));
    Rfan_tip2 = Rfan_tip1 * 0.98;
     
    
    %% FAN/IPC transition ducts
    Rtrans_duct_hub_1 = Rfan_hub2;
    AR_trans_duct = 0.4;

    v_2 = Rfan_hub2/Rfan_tip2;
    A2 = pi*Rfan_tip2^2 * (1 - v_2^2);
    A_duct_entry = A2/ (BPR + 1);
    r_spiltter_lip = sqrt((A_duct_entry/pi) + Rtrans_duct_hub_1^2);
    l_ax_trans_duct = (r_spiltter_lip - Rtrans_duct_hub_1)/AR_trans_duct;
    

    %% IPC
    M_ax_IPC_in = 0.539;
    M_ax_IPC_out = 0.341;
    M_ax_IPC_relative = 1.55;
    v_IPC_in = 0.630;
    v_IPC_out = 0.819;


    %IPC Rhub and Rtip 1
    IPCArea_1 = calculateFanArea(M_ax_IPC_in ,gamma_air,T_02,P_02,m_h);
    RIPC_tip1 = sqrt(IPCArea_1/(pi*(1-v_IPC_in^2)));
    RIPC_hub1 = v_IPC_in * RIPC_tip1;
 
        %The exit hub and tip radius are set equal to the entry hub and tip radius of the IPC.
        Rtrans_duct_hub_2 = RIPC_hub1;
        Rtrans_duct_tip_2 = RIPC_tip1;

     %IPC Rhub and Rtip 2
    IPCArea_2 = calculateFanArea(M_ax_IPC_out ,gamma_air,T_03,P_03,m_h);
    RIPC_tip2 = sqrt(IPCArea_2/(pi*(1-v_IPC_out^2)));
    RIPC_hub2 = v_IPC_out * RIPC_tip2;

    %Average stage loadings and stage numbers
    % Constants and parameters
      EIS_IPC = -8.968 + (0.004877 * 2020); 
      R_mid_1 = (RIPC_tip1 + RIPC_hub1)/2;
      R_mid_2 = (RIPC_tip2 + RIPC_hub2)/2;
      N = 2;
      T_2 = T_02/(1+M_ax_IPC_in^2*(gamma_air-1)/2);

      delta_H = cp_air * (T_03-T_02);
      C_a_IPC_tip_1 = M_ax_IPC_relative * sqrt(gamma_air*R*T_2);
      C_a_IPC_tip_accual = M_ax_IPC_in * sqrt(gamma_air*R*T_2);
      U_IPC_tip = sqrt(C_a_IPC_tip_1^2 - C_a_IPC_tip_accual^2); %Velocity triangle
      Omega_IPC_rotor = U_IPC_tip / RIPC_tip1;
      U_mid_1 = Omega_IPC_rotor * R_mid_1;
      U_mid_last_stage = Omega_IPC_rotor * R_mid_2;
      U_mid = U_mid_1 + (U_mid_last_stage - U_mid_1) * (0:N-1) / (N-1);
      avg_stage_loading = 2*delta_H/(sum(U_mid.^2));

      while avg_stage_loading >= EIS_IPC
      N = N + 1;
      U_mid = U_mid_1 + (U_mid_last_stage - U_mid_1) * (0:N-1) / (N-1);
      avg_stage_loading = 2*delta_H/(sum(U_mid.^2));
      end

     %Component length – rapid procedure
     c_IPC = 0.3;
     AR_IPC_entry = 36.20 - (0.01694 * EIS);
     AR_IPC_exit = 35.47 - (0.01694 * EIS);
     h1_IPC = RIPC_tip1 - RIPC_hub1;
     h2_IPC = RIPC_tip2 - RIPC_hub2;
     h_hat_IPC = sqrt(h1_IPC*h2_IPC);
     l_ax_IPC = 2 * N * (h_hat_IPC/((AR_IPC_exit + AR_IPC_entry)/2)) * (1 + c_IPC);



     %% HPC
     M_ax_HPC_in = 0.482;
     M_ax_HPC_out = 0.263;
     M_ax_HPC_relative = 1.3;
     v_HPC_in = 0.5613/(1.487-exp(-0.04286*(HPC+0.5718)));
     v_HPC_out = 0.908;
        
     %HPC Rhub and Rtip 1  
     HPCArea_1 = calculateFanArea(M_ax_HPC_in ,gamma_air,T_03,P_03,m_h);
     RHPC_tip1 = sqrt(HPCArea_1/(pi*(1-v_HPC_in^2)));
     RHPC_hub1 = v_HPC_in * RHPC_tip1;

     %HPC Rhub and Rtip 2
     HPCArea_2 = calculateFanArea(M_ax_HPC_out ,gamma_air,T_04,P_04,m_h_corr);
     RHPC_tip2 = sqrt(HPCArea_2/(pi*(1-v_HPC_out^2)));
     RHPC_hub2 = v_HPC_out * RHPC_tip2;


     %Average stage loadings and stage numbers
     % Constants and parameters
      EIS_HPC = -5.736 + 0.00323 * EIS;
      R_mid_HPC_1 = (RHPC_tip1 + RHPC_hub1)/2;
      R_mid_HPC_2 = (RHPC_tip2 + RHPC_hub2)/2;
      N_HPC = 2;
      delta_H_HPC = cp_air * (T_04-T_03);
      T_3 = T_03/(1+M_ax_HPC_in^2*(gamma_air-1)/2);
      C_a_HPC_tip_1 = M_ax_HPC_relative * sqrt(gamma_air*R*T_3);
      C_a_HPC_tip_accual = M_ax_HPC_in * sqrt(gamma_air*R*T_3);
      U_HPC_tip = sqrt(C_a_HPC_tip_1^2 - C_a_HPC_tip_accual^2); %Velocity triangle
      Omega_HPC_rotor = U_HPC_tip / RHPC_tip1;
      U_mid_HPC_1 = Omega_HPC_rotor * R_mid_HPC_1;
      U_mid_HPC_last_stage = Omega_HPC_rotor * R_mid_HPC_2;
      U_mid_HPC = U_mid_HPC_1 + (U_mid_HPC_last_stage - U_mid_HPC_1) * (0:N_HPC-1) / (N_HPC-1);
      avg_stage_loading_HPC = 2*delta_H_HPC/(sum(U_mid_HPC.^2));
      
     
      while avg_stage_loading_HPC > EIS_HPC
      N_HPC = N_HPC + 1;
      U_mid_HPC = U_mid_HPC_1 + (U_mid_HPC_last_stage - U_mid_HPC_1) * (0:N_HPC-1) / (N_HPC-1);
      avg_stage_loading_HPC = 2*delta_H_HPC/(sum(U_mid_HPC.^2));
      end
    
      %Component length – rapid procedure
      c_HPC = 0.3;
      AR_HPC_entry = 31.4 - (0.0147 * EIS);
      AR_HPC_exit = 30.7 - (0.0147 * EIS);
      h1_HPC = RHPC_tip1 - RHPC_hub1;
      h2_HPC = RHPC_tip2 - RHPC_hub2;
      h_hat_HPC = sqrt(h1_HPC*h2_HPC);
      l_ax_HPC = 2 * N_HPC * (h_hat_HPC/((AR_HPC_exit + AR_HPC_entry)/2)) * (1 + c_HPC);


       %% IPC/HPC transsition duct
     Rtrans_duct_IPCHPC_tip_1 = RIPC_tip2;
     Rtrans_duct_IPCHPC_hub_1 = RIPC_hub2;
     Rtrans_duct_IPCHPC_hub_2 = RHPC_tip1;
     Rtrans_duct_IPCHPC_tip_2 = RHPC_hub1;
     l_ax_trans_duct_IPCHPC = (Rtrans_duct_IPCHPC_tip_1 - Rtrans_duct_IPCHPC_hub_1)/AR_trans_duct;
      

      %% Combustor
       %Assumed combustor average Mach number is 0.06
       % , beta is 10 deg and time is 6 ms
       M_avg_combustor = 0.06;
       beta_combustor = 10;
       t_combustor = 6 * 10^-3;

       T_combustor_avg = (T_04 + T_05)/2;
       P_combustor_avg = (P_04 + P_05)/2;

       COMBUSTORArea = calculateFanArea(M_avg_combustor,gamma_gas,T_combustor_avg,P_combustor_avg,(m_h_corr + m_fuel));

       T_4 = T_combustor_avg/(1+M_avg_combustor^2*(gamma_air-1)/2);

      

       V_combustor = M_avg_combustor * sqrt(gamma_gas*R*T_4); 
       l_ax_combustor = t_combustor * V_combustor;
       dy_combustor = l_ax_combustor * tand(beta_combustor);
      

       %% HPT
       EIS_HPT_IPT = 3.247;
       c_HPT = 0.2;
       AR_HPT = 29.233 - (0.0140*EIS);

       M_ax_HPT_in = 0.150;
       M_ax_HPT_out = 0.331 + (0.0610*(P_05/P_06));

       R_mid_HPT_1 = R_mid_HPC_2 + dy_combustor;
       R_mid_HPT_2 = R_mid_HPT_1;

       HPTArea_1 = calculateFanArea(M_ax_HPT_in ,gamma_gas,T_05,P_05,(m_h_corr+m_fuel));
       HPTArea_2 = calculateFanArea(M_ax_HPT_out ,gamma_gas,T_06,P_06,(m_h+m_fuel));

        h1_HPT = HPTArea_1/(R_mid_HPT_1*2*pi); %annulus area formula
        RHPT_tip1 =  R_mid_HPT_1  + (h1_HPT/2);
        RHPT_hub1 =  R_mid_HPT_1  - (h1_HPT/2);

        h2_HPT = HPTArea_2/(R_mid_HPT_2*2*pi); %annulus area formula
        RHPT_tip2 =  R_mid_HPT_2  + (h2_HPT/2);
        RHPT_hub2 =  R_mid_HPT_2  - (h2_HPT/2);

        h_hat_HPT = sqrt(h1_HPT*h2_HPT);

        l_ax_HPT = 2 * (h_hat_HPT/AR_HPT) * (1 + c_HPT);
      
        AN2_HPT = HPTArea_1*(Omega_HPC_rotor/(pi*2))^2; 

        AN2_HPT_EIS = (-158.99 + 2020 * 0.0830) * 10^3;

        avg_stage_loading_HPT = 2*(delta_H_HPC/Shaft_mec_eff)/ (Omega_HPC_rotor * R_mid_HPT_1)^2;

        if AN2_HPT <= AN2_HPT_EIS
            disp('HPT STRESS OK!')
        end

        if avg_stage_loading_HPT <= EIS_HPT_IPT
            disp('HPT STAGE LOAD OK!')
        end

        




        %% IPT
        AR_IPT = 30.492 - (0.0140 * 2020);
        M_ax_IPT_in = 0.314;
        M_ax_IPT_out = 0.554;
        
        IPTArea_1 = calculateFanArea(M_ax_IPT_in ,gamma_gas,T_06,P_06,(m_h+m_fuel));
        RIPT_tip1 = 1.3 * RHPT_tip2; %From duct explantion 
        RIPT_hub1 = sqrt(RIPT_tip1^2 - (IPTArea_1/pi)); %From pi*r_t^2-pi*r_h^2 == HPTArea_1;
        R_mid_IPT_1 = (RIPT_tip1 + RIPT_hub1) / 2;
            

        IPTArea_2 = calculateFanArea(M_ax_IPT_out ,gamma_gas,T_07,P_07,(m_h+m_fuel));
        R_mid_IPT_2 = R_mid_IPT_1; %Assume same mid radius
        h2_IPT = IPTArea_2/(R_mid_IPT_2*2*pi); %annulus area formula
        RIPT_tip2 =  R_mid_IPT_2  + (h2_IPT/2);
        RIPT_hub2 =  R_mid_IPT_2  - (h2_IPT/2);

        avg_stage_loading_IPT = 2*(delta_H/Shaft_mec_eff)/ (Omega_IPC_rotor * R_mid_IPT_1)^2;

        h1_IPT = RIPT_tip1 - RIPT_hub1;
        h2_IPT = RIPT_tip2 - RIPT_hub2;
        h_hat_IPT = sqrt(h1_IPT*h2_IPT);

        l_ax_IPT = 2 * (h_hat_IPT/AR_IPT) * (1 + c_HPT);
      
        AN2_IPT = IPTArea_1*(Omega_IPC_rotor/(pi*2))^2; 

        AN2_IPT_EIS = (-133.97 + (2020 * 0.07099)) * 10^3;

        if AN2_IPT <= AN2_IPT_EIS
            disp('IPT STRESS OK!')
        else
            disp('IPT STRESS NOT OK!')
        end

        if avg_stage_loading_IPT <= EIS_HPT_IPT
            disp('IPT STAGE LOAD OK!')
        else 
             disp('IPT STAGE LOAD NOT OK!')
        end


       
        
        %% HPT/IPT Transsition ducts
        Rtrans_duct_HPTIPT_tip_1 = RHPT_tip2;
        Rtrans_duct_HPTIPT_hub_1 = RHPT_hub2;
        Rtrans_duct_HPTIPT_hub_2 = RIPT_hub1;
        Rtrans_duct_HPTIPT_tip_2 = RIPT_tip1;
        l_ax_trans_duct_HPTIPT = (Rtrans_duct_HPTIPT_tip_1 - Rtrans_duct_HPTIPT_hub_1)/AR_trans_duct;


        %% LPT
        EIS_LPT = -39.26 + (0.02185 * 2020);
        c_LPT = 0.4;
        AR_LPT_entry = 31.40 - (0.0146 * 2020);
        AR_LPT_exit = -35.60 + (0.0209 * 2020);

        M_ax_LPT_in = 0.368;
        M_ax_LPT_out = 0.322;

        LPTArea_1 = calculateFanArea(M_ax_LPT_in ,gamma_gas,T_07,P_07,(m_h+m_fuel));
        RLPT_tip1 = 1.15 * RIPT_tip2; %From duct explantion 
        RLPT_hub1 = sqrt(RLPT_tip1^2 - (LPTArea_1/pi)); %From pi*r_t^2-pi*r_h^2 == HPTArea_1;
        R_mid_LPT_1 = (RLPT_tip1 + RLPT_hub1) / 2;

        LPTArea_2 = calculateFanArea(M_ax_LPT_out ,gamma_gas,T_08,P_08,(m_h+m_fuel));
        RLPT_tip2 = Rfan_tip1 * (1.25/(1.976-exp(-0.2503*(BPR+0.6410))));
        RLPT_hub2 = sqrt(RLPT_tip2^2 - (LPTArea_2/pi));
        R_mid_LPT_2 = (RLPT_tip2 + RLPT_hub2) /2;

        %Average stage loadings and stage numbers
        N_LPT = 1;
        delta_H_LPT = (cp_air * (T_07-T_08))/Shaft_mec_eff;
        U_fan_tip = (-(59.74*FPR) + (88.07*(FPR^2)) - (25.93*(FPR^3)))*sqrt(T_01);
          Omega_LPT_rotor = U_fan_tip / Rfan_tip1;
          U_mid_1_LPT = Omega_LPT_rotor * R_mid_LPT_1;
          U_mid_last_stage_LPT = Omega_LPT_rotor * R_mid_LPT_2;
          U_mid_LPT = U_mid_1_LPT^2;
          avg_stage_loading_LPT = 2*delta_H_LPT/(sum(U_mid_LPT));

      while avg_stage_loading_LPT >= EIS_LPT
          N_LPT = N_LPT + 1;
          U_mid_LPT = U_mid_1_LPT + (U_mid_last_stage_LPT - U_mid_1_LPT) * ((0:N_LPT-1) / (N_LPT-1));
          U_mid_LPT = U_mid_LPT.^2;
          avg_stage_loading_LPT = 2*delta_H_LPT/(sum(U_mid_LPT));
      end
        
     %Component length – rapid procedure
     h1_LPT = RLPT_tip1 - RLPT_hub1;
     h2_LPT = RLPT_tip2 - RLPT_hub2;
     h_hat_LPT = sqrt(h1_LPT*h2_LPT);
     l_ax_LPT = 2 * N_LPT * (h_hat_LPT/((AR_LPT_exit + AR_LPT_entry)/2)) * (1 + c_LPT);


     AN2_LPT = LPTArea_2*(Omega_LPT_rotor/(pi*2))^2;


     %% IPT/LPT Transsition ducts
        Rtrans_duct_IPTLPT_tip_1 = RIPT_tip2;
        Rtrans_duct_IPTLPT_hub_1 = RIPT_hub2;
        Rtrans_duct_IPTLPT_hub_2 = RLPT_hub1;
        Rtrans_duct_IPTLPT_tip_2 = RLPT_tip1;
        l_ax_trans_duct_IPTLPT = (Rtrans_duct_IPTLPT_tip_1 - Rtrans_duct_IPTLPT_hub_1)/AR_trans_duct;

     y_values = [Rfan_hub1, Rfan_tip1, Rfan_tip2,Rfan_hub2;   % y-values for component 1
                 Rtrans_duct_hub_1, r_spiltter_lip, Rtrans_duct_tip_2,Rtrans_duct_hub_2;
                 RIPC_hub1,RIPC_tip1,RIPC_tip2,RIPC_hub2;
                 Rtrans_duct_IPCHPC_hub_1,Rtrans_duct_IPCHPC_tip_1,Rtrans_duct_IPCHPC_hub_2,Rtrans_duct_IPCHPC_tip_2;
                 RHPC_hub1,RHPC_tip1,RHPC_tip2,RHPC_hub2;
                 RHPC_hub2,RHPC_tip2,RHPT_tip1,RHPT_hub1;
                 RHPT_hub1,RHPT_tip1, RHPT_tip2, RHPT_hub2;
                 Rtrans_duct_HPTIPT_hub_1,Rtrans_duct_HPTIPT_tip_1, Rtrans_duct_HPTIPT_tip_2,Rtrans_duct_HPTIPT_hub_2;
                 RIPT_hub1,RIPT_tip1, RIPT_tip2, RIPT_hub2;
                 Rtrans_duct_IPTLPT_hub_1,Rtrans_duct_IPTLPT_tip_1,Rtrans_duct_IPTLPT_tip_2,Rtrans_duct_IPTLPT_hub_2;
                 RLPT_hub1, RLPT_tip1, RLPT_tip2, RLPT_hub2];  
    

    % Define x-distances for each component
    x_distances = [0, l_ax_fan;   % x-distances for component 1
                    l_ax_fan,l_ax_fan+l_ax_trans_duct;
                    l_ax_fan+l_ax_trans_duct,l_ax_fan+l_ax_trans_duct + l_ax_IPC;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT + l_ax_IPT;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT + l_ax_IPT, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT + l_ax_IPT + l_ax_trans_duct_IPTLPT;
                    l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT + l_ax_IPT + l_ax_trans_duct_IPTLPT, l_ax_fan+l_ax_trans_duct + l_ax_IPC + l_ax_trans_duct_IPCHPC + l_ax_HPC + l_ax_combustor + l_ax_HPT + l_ax_trans_duct_HPTIPT + l_ax_IPT + l_ax_trans_duct_IPTLPT + l_ax_LPT];  % x-distances for component 2 (example, add more rows for additional components) 

    %sketch(y_values,x_distances);

    rearranged_values = y_values(:, [1, 2, 4, 3]);

    differences = x_distances(:, 2) - x_distances(:, 1);






%% D3


 %Station 01
    M_ax_01 = M_ax_HPC_in;
    Area_rotor_in = HPCArea_1;
    tin = T_03/(1+M_ax_01^2*(gamma_air-1)/2);
    vin = sqrt(gamma_air*R*tin);
    cin = vin*M_ax_01;

    Tin = T_03;
    Pin = P_03;
    mf = m_h;
    h1 = h1_HPC;


  %Station 2
   %Calculating stagnation temperature and pressure
   T02 = ((avg_stage_loading_HPC*U_mid_HPC(1)^2/2) + (cp_air*T_03))/cp_air;
   P02 = ((T02/T_03)^((HPC_poly_eff*gamma_air)/(gamma_air-1)))* P_03;
   

   %Calculating the area at the end of the stage
   mach_numbers = linspace(M_ax_HPC_in, M_ax_HPC_out, N_HPC*2);
   M_ax_2 = mach_numbers(2);
   Area_stator_out = calculateFanArea(M_ax_2,gamma_air,T02,P02,m_h);

   tout = P02/(1+M_ax_2^2*(gamma_air-1)/2);
   vout = sqrt(gamma_air*R*tout);
   cout = vout*M_ax_2;

   %Calculating the height at the end of the stage
   mid_radi = linspace(R_mid_HPC_1, R_mid_HPC_2, N_HPC*2);
   mid_radi_end = mid_radi(2);
   h2 = Area_stator_out/(mid_radi_end*2*pi); 
        rstat_tte =  mid_radi_end  + (h2/2);
        rstat_hte =  mid_radi_end  - (h2/2);
    

   %Station 12
    h_hat = sqrt(h1*h2);
    AR = linspace(AR_HPC_entry, AR_HPC_exit, N_HPC*2); 

    AR_rotor = AR(1); 
    AR_stator = AR(2);
        
    l_rotor = h_hat/AR_rotor;
    l_stator = h_hat/AR_stator;
    l_space_rtst = l_rotor * 0.3;
    l_igv = 0.8 * l_rotor;
    l_igv_space = 0.3*l_igv;

    x2 = (l_igv_space+l_rotor+l_space_rtst+l_stator);
    rrot_tle = radii_interpol(RHPC_tip1,RHPC_tip1,rstat_tte,0,x2,(l_igv_space));
    rrot_tte = radii_interpol(RHPC_tip1,RHPC_tip1,rstat_tte,0,x2,(l_igv_space+l_rotor));
    rstat_tle = radii_interpol(RHPC_tip1,RHPC_tip1,rstat_tte,0,x2,(l_igv_space+l_rotor+l_space_rtst));

    rrot_hle = radii_interpol(RHPC_hub1,RHPC_hub1,rstat_hte,0,x2,(l_igv_space));
    rrot_hte = radii_interpol(RHPC_hub1,RHPC_hub1,rstat_hte,0,x2,(l_igv_space+l_rotor));
    rstat_hle = radii_interpol(RHPC_hub1,RHPC_hub1,rstat_hte,0,x2,(l_igv_space+l_rotor+l_space_rtst));
    
 %Stattion 0
    %Radial
    r_igv_tle = RHPC_tip1;
    r_igv_hle = RHPC_hub1;
    r_igv_tte = r_igv_tle;
    r_igv_hte = r_igv_hle;

 %Axial coordinates
 x_igv_le = 0;
 x_igv_te = l_igv;
 x_rot_le = x_igv_te + l_igv_space;
 x_rot_te = x_rot_le + l_rotor;
 x_stat_le = x_rot_te + l_space_rtst;
 x_stat_te = x_stat_le + l_stator;
     
%saveCompressorForSmooth((Omega_HPC_rotor/(2*pi)), m_h, tin, Pin, P02, x_igv_le, x_igv_te, x_rot_le,x_rot_te, x_stat_le, x_stat_te, r_igv_tle, r_igv_tte, rrot_tle, rrot_tte, rstat_tle,rstat_tte, r_igv_hle, r_igv_hte, rrot_hle, rrot_hte, rstat_hle, rstat_hte);




    




    
      
