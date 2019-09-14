%=========================================================================
% FILENAME:     RocketSizing.m
% WRITTEN BY:   Dakota Jandek
%=========================================================================

clc;
clear;
format short g;
%% Init

dt = .0001;                     %[s]        Time Step

g = 32.2;                       %[ft/s/s]   Newton Gravitational Constant
Isp = 165;                      %[s]        Specific Impulse
m_prop = 27;                    %[lbm]      Propellant Mass
m_dry = 125;                    %[lbm]      Dry Mass
mdot_inital = 9.39;             %[lbm/s]    Inital Mass Flow Rate
mdot_average = 4.21;            %[lbm/s]    Average Mass Flow Rate
tb = m_prop/mdot_average;       %[s]        Burn Time
rocket_radius = .25;            %[ft]       Radius of the Rocket
Cd = .5;                        %[]         Drag Coefficient

mdot = linspace(mdot_inital, mdot_inital-(2*(mdot_inital-mdot_average)), ceil(tb));     %[lbm/s] Linear Fit Mass Flow Rate
T = Isp*mdot;                   %[lbm]      Thrust Curve
I = Isp*mdot_average*tb;        %[lbm]      Total Impulse   

launchRailLength = 40;          %[ft]       Length of the Launch Rail

%Iteration variables
m0 = m_dry + m_prop;                       
t = 0;
v = 0;
h = 0;
v_t(1) = 0;
h_t(1) = 0;
t_t(1) = 0;
D_t(1) = 0;
gt_t(1) = 0;
t_t(1) = 0;
m1 = m0;
i = 1;
bool_maxVel = 0;
bool_rail = 0;

v_bo = 0;                       %[ft/s]     Velocity of Burn Out
velocityAtLaunchRail = 0;       %[ft/s]     Velocity at the End of the Launch Rail
h_bo = 0;                       %[ft]       Height of Burn Out
t_apogee = 0;                   %[s]        Time of Apogee
%% Iteration

while v >= 0
    if(bool_rail == 0 && h >= launchRailLength)
        velocityAtLaunchRail = v;
        bool_rail = 1;
    end
    if(m1 > m_dry) %There is still propellent to be burned
        mdot_current = mdot(floor(t) + 1);
        m2 = m1;
        m1 = m1 - mdot_current*dt;
        v = v + g*Isp*log(m2/m1) - g*dt - Cd*.5*getRho(h)*v^2*(pi*rocket_radius^2)/(m2)*dt;       
    else %All propellent has been burned, maximum velocity has been reached
        if(bool_maxVel == 0)
            v_bo = v;
            bool_maxVel = 1;
            h_bo = h;
        end
        v = v - g*dt - Cd*.5*getRho(h)*v^2*(pi*rocket_radius^2)/(m2)*dt;
    end
    h = h + v*dt;
    i = i+1;
    t = t + dt;
    t_t(i) = t;
    v_t(i) = v;
    h_t(i) = h;
    t_t(i) = t;
    D_t(i) = Cd*.5*getRho(h)*v^2*(pi*rocket_radius^2)/((m2+m1)/2)*dt;
    gt_t(i) = g*dt;
end

t_apogee = t_t(length(t_t));
%% Plotting
figure(1)

subplot(4,1,1)
plot(t_t, h_t);
xlabel('Time (s)');
ylabel('Height (ft)');
title('Height(t)');

subplot(4,1,2)
plot(t_t, v_t);
xlabel('Time (s)');
ylabel('Velocity (ft/s)');
title('Velocity(t)');

subplot(4,1,3)
plot(t_t, gt_t/dt);
xlabel('Time (s)');
ylabel('Gravity Loss (ft)');
title('Gravity Loss(t)');

subplot(4,1,4)
plot(t_t, D_t/dt);
xlabel('Time (s)');
ylabel('Drag Loss(ft)');
title('Drag Loss(t)');

%% Command Window Output
fprintf('Dry Mass: %35.0f lbm\n',m_dry);
fprintf('Isp: %40.0f s\n',Isp);
fprintf('Average Mass Flow: %26.2f lbm/s\n',mdot_average);
fprintf('Gravitational Constant: %21.3f ft/s/s\n',g);
fprintf('Assumed mass of propellent: %17.3f lbm\n',m_prop);
fprintf('Final Height: %31.3f ft\n',h_t(i));
fprintf('Burn Time: %34.3f s\n',tb);
fprintf('Time to Apogee: %29.3f s\n',t_apogee);
fprintf('Velocity at end of Launch Rail: %13.3f ft/s\n',velocityAtLaunchRail);
fprintf('Maximum Velocity Reached: %19.3f ft/s\n',v_bo);
fprintf('Total Loss to Drag: %25.3f ft\n',sum(D_t)*t_apogee);
fprintf('Total Loss to Gravity: %22.3f ft\n\n',sum(gt_t)*t_apogee);
fprintf('Average Thrust: %29.3f lbm\n',Isp*mdot_average);
fprintf('Total Impulse: %30.3f lbm\n',I);





%% Density of Air as a Function of Altitude
function rho = getRho(h)
    switch h
        case h + 4595 < 5000
            rho = .0659;
        case h  + 4595 < 10000
            rho = .0565;
        case h + 4595 < 15000
            rho = .0481;
        otherwise 
            rho = .0408;
    end
end




