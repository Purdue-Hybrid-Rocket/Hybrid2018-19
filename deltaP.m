%% Change the following
% line 10: mass flow rate 
% line 12: pipe diameter 
% line 15: n2o density in the pipe
% line 16: get number from Paul 
% line 22: friction factor (depends on surface roughness and inlet/outlet shape)
% line 27: valve cv
%% Pressure drop calculation (Needs Improvement)

mdot = 4.55*0.453592; %kg/s
alpha = 1;
pipe_D = .402; %in  
A_pipe = (pipe_D/2 * .0254)^2 * pi; %m^2
rho_pipe = 52.92*16.0185; %kg/m3
V_A = mdot ./ (rho_pipe*A_pipe); %m/s
g = 8.9*9.8; %8.9g, m/s
P_A = 1000*6894.76; %Pa

%A to B (Tank to Main Valve Entrance)
L_AB = 7; %in
delZ_AB = 7; %in
f = .03;
h_LT_AB = f * L_AB/pipe_D * V_A.^2/(2*g); %m
P_B = P_A + (rho_pipe .* g .* (delZ_AB .* .0254 - h_LT_AB))/1000000; %MPa

%B to C (Main Valve Entrance to Main Valve Exit)
cv = 12;
q = mdot./rho_pipe*15850.372483753; %GPM %.5 * rho_pipe .* V_A.^2; %Pa
SG = rho_pipe/997;
delP_BC = q.^2./cv^2 .* SG; %psi
P_C = P_B - delP_BC/145.038; %MPa

%C to D (Main Valve Exit to Injector Plate)
L_CD = 6.5; %in
delZ_CD = 6.5; %in
f = .03;
h_LT_CD = f * L_CD/pipe_D * V_A.^2/(2*g); %m
P_D = P_C + (rho_pipe .* g .* (delZ_CD .* .0254 - h_LT_CD))/1000000; %MPa

%pressure drops (psi)
PA = (P_A-P_B)*145.038
PB = (P_B-P_C)*145.038
PC = (P_C-P_D)*145.038
delta_P = PA + PB + PC 

%velocity
velocity = ((4.55*4) / (.402^2 * pi * 52.92 * .0005787)) * 0.0833333 %ft/s