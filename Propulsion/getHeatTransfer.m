function [max_time] = getHeatTransfer(Ti_list,chamber_length,h_list,k_list,c_list,t_list,rho_list,T_max,time_step,total_time)

%%This function simulates transient heat flow of a fluid within a combustion chamber
%%with various layers of solid material. Internal fluid is assumed as infinite
%%reservoir (no change in temperature) 
%Input descriptions:
%Ti_list = Initial temperature list, including ambient air [K]
%chamber_length = combustion chamber section length [in]
%h_list = fluid convection coefs [W/(cm^2-K)]
%k_list = material conductive coefs [W/(cm-K)]
%c_list = heat capacity coefs [J/(g-K)] (first and last entry zeros)
%t_list = layer thicknesses [in] (first entry is chamber radius, include ambient air)
%rho_list = density list [g/cm^3]

%Example Input:
%[max_time] = getHeatTransfer([3300 295 295 295],2,[1.0758 0 0 0.00026],[0 0.00293 2.36 0],[1 1.3 0.921 1.005],[2.625 0.25 0.125 2],[1 1.8 2.7 0.001225],400,0.01,9)

%% INITIALIZATION %%
total_i = floor(total_time / time_step) + 1; %total number of iterations
total_layers = length(t_list); %total number of layers including fluid
T_hist = zeros(total_i,total_layers); %temperature history of each layer
T_hist(:,1) = Ti_list(1);
T_hist(:,4) = Ti_list(4);
T_hist(1,:) = Ti_list;
time_hist = linspace(0,total_time,total_i);
t_list = t_list * 2.54; %convert from [in] to [cm]
chamber_length = chamber_length * 2.54; %convert from [in] to [cm]
r_list = t_list; %radius list [cm]
for L = 1:(total_layers-1)
r_list(L+1) = r_list(L) + t_list(L+1);
end
A_list = 2 * pi * r_list * chamber_length; %wall surface area [cm^2]
m_list(1) = pi * r_list(1) ^ 2 * chamber_length * rho_list(1); %mass of each layer [g]
for L = 2:total_layers
m_list(L) = pi * (r_list(L) ^ 2 - r_list(L-1) ^ 2) * chamber_length * rho_list(L);
end
T_max = zeros(total_i,1) + T_max;
%% MAIN LOOP %%
for i = 1:(total_i - 1)
    for L = 1:(total_layers - 1)
        if h_list(L) > 0 && k_list(L+1) > 0 %fluid-solid boundary
            q_dot = h_list(L) * (T_hist(i+1,L) - T_hist(i,L+1)); %heat flux [W/cm^2]
            delta_Q = A_list(L) * q_dot * time_step;
            delta_T = delta_Q / (m_list(L+1) * c_list(L+1)); %solid layer temp change
            T_hist(i+1,L+1) = T_hist(i,L+1) + delta_T;

            if L > 1 %conserve energy after first layer
                delta_T = -delta_Q / (m_list(L) * c_list(L)); %fluid layer temp change
                T_hist(i+1,L) = T_hist(i+1,L) + delta_T;
            end

        elseif k_list(L) > 0 && h_list(L+1) > 0 %solid-fluid boundary
            %%q_dot = 0;
            q_dot = h_list(L+1) * (T_hist(i+1,L) - T_hist(i,L+1)); %heat flux [W/cm^2]
            delta_Q = A_list(L) * q_dot * time_step;
            %%delta_T = delta_Q / (m_list(L+1) * c_list(L+1)); %fluid layer temp change
            %%T_hist(i+1,L+1) = T_hist(i,L+1) + delta_T;

            delta_T = -delta_Q / (m_list(L) * c_list(L)); %solid layer temp change
            T_hist(i+1,L) = T_hist(i+1,L) + delta_T;

        elseif k_list(L) > 0 && k_list(L+1) > 0 %solid-solid boundary
            q_dot = k_list(L) / t_list(L) * (T_hist(i+1,L) - T_hist(i,L+1)); %heat flux [W/cm^2]
            delta_Q = A_list(L) * q_dot * time_step;
            delta_T = delta_Q / (m_list(L+1) * c_list(L+1)); %second solid layer temp change
            T_hist(i+1,L+1) = T_hist(i,L+1) + delta_T;

            delta_T = -delta_Q / (m_list(L) * c_list(L)); %first solid layer temp change
            T_hist(i+1,L) = T_hist(i+1,L) + delta_T;

        else
            sprintf('Error: invalid boundary condition')
            return
        end
    end
end
%% OUTPUT %%
for i = 1:total_i
    if T_hist(i,3) > T_max
    max_time = i / total_i * total_time;
    break
    else
    max_time = 0;
    end
end
T_hist = (T_hist - 273.15) * 9 / 5 + 32;
plot(time_hist,T_hist(:,1),'y-',time_hist,T_hist(:,2),'r-',time_hist,T_hist(:,3),'k-',time_hist,T_hist(:,4),'b-')
title('0.375" Thick Insulation');
xlabel('Exposure Time [s]');
ylabel('Temperature [F]');
legend('Chamber','Inner Insulation','Inner Aluminum','Ambient Air');
set(gcf,'color','w');
grid on;
