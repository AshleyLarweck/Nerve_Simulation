%===============================  NERVE ACTION POTENTIAL SIMULATION ========================================
%============================================================================================================
close all;clear all;clc;
%============================================================================================================
%                                         Given Parameters for simulation
%============================================================================================================
E_rest=-68;
E_K=-82;
E_Na=45;
C=1;
G_K=36;
G_Na=120;
n_init=0.2508;
m_init=0.5508
h_init=0.1534;
%==========================================================================================================
%                                          Determining the Initial values
%                                             Set the initial states
%===========================================================================================================
V =0; % baseline voltage
alpha_n=0.01*(10-V)/(exp((10-V)/10)-1);
beta_n=0.125*exp(-V/80);
alpha_m=0.1*(25-V)/(exp((25-V)/10)-1);
beta_m=4*exp(-V/18);
alpha_h=0.07*exp(-V/20);
beta_h=1/(exp((30-V)/10)+1);
n(1)=alpha_n/(alpha_n+beta_n);
m(1)=alpha_m/(alpha_m+beta_m);
h(1)=alpha_h/(alpha_h+beta_h);
%=========================================================================================================
%                                  Defining the time vector for simulation
%===========================================================================================================
Total_time=50;  %50 ms
deltaT=0.01;
t=0:deltaT:Total_time;
%============================================================================================================
%                                               Defining the Current
%===============================================================================================================
I=zeros(1,numel(t));
I(1,1:80)=75; % % micro amp current

%===========================================================================================================
%                            Computing coefficients, currents, and derivatives at each time step
%===========================================================================================================
for i=1:numel(t)-1

    %---calculate the coefficients---%

    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = .01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) );
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = .1*( (25-V(i)) / (exp((25-V(i))/10)-1) );
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);


    %---calculate the currents---%

    I_Na = (m(i)^3) *G_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K = (n(i)^4) * G_K * (V(i)-E_K); %Equations 4 and 6
    I_ion = I(i) - I_K - I_Na ;

    %---calculate the derivatives using Euler first order approximation---%
    V(i+1) = V(i) + deltaT*I_ion/C;
    n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i)); %Equation 7
    m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i)); %Equation 15
    h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i)); %Equation 16

end

V = V-70; %Set resting potential to -70mv


%============================================================================================================
                                                 % Plot the Voltage
%============================================================================================================
plot(t,V,'LineWidth',3)
hold on
legend({'Voltage'})
ylabel('Voltage (mV)')
xlabel('time (ms)')
title('Voltage over time in simulated neuron');
%============================================================================================================
                                                 %Plot the Conductance
%============================================================================================================
figure
p1 = plot(t,G_K*n.^4,'m','LineWidth',2);
hold on
p2 = plot(t,G_Na*(m.^3).*h,'g','LineWidth',2);
legend([p1, p2], 'Conductance for Potassium', 'Conductance for Sodium')
ylabel('Conductance')
xlabel('time (ms)')
title('Conductance for Potassium and Sodium Ions in Simulated Neuron')




function HodgkinHuxley
Vrest = -70; % mV 
dt = 0.01; % ms
totalTime = 50; % ms
C = 1.0; % microF/cm^2
% constants
E_Na = 115 + Vrest; % mV
E_K = -6 + Vrest; %mV
E_Leak = 10.6 + Vrest; % mV
g_Na = 120; % mS/cm^2
g_K = 36; % mS/cm^2
g_Leak = 0.3; % mS/cm^2
% matrix of timesteps
t = [0:dt:totalTime];
% Current input −− change this to see how different inputs affect the neuron
I_current = zeros(1,length(t));
I_current(10/dt:11) = 15; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time


% initializing values
V(1) = Vrest; % membrane potential is starting at rest
% calculation of alpha and beta values
[alphaM, betaM] = m_equations(V(1), Vrest);
[alphaN, betaN] = n_equations(V(1), Vrest);
[alphaH, betaH] = h_equations(V(1), Vrest);
% initializing gating variables to the asymptotic values when membrane potential
% is set to the membrane resting value based on equation 13
m(1) = (alphaM / (alphaM + betaM));
n(1) = (alphaN / (alphaN + betaN));
h(1) = (alphaH / (alphaH + betaH));
% repeat for time determined in totalTime , by each dt
for i = 1:length(t)
% calculate new alpha and beta based on last known membrane potenatial
[alphaN, betaN] = n_equations(V(i), Vrest);
[alphaM, betaM] = m_equations(V(i), Vrest);
[alphaH, betaH] = h_equations(V(i), Vrest);
% conductance variables − computed separately to show how this
% changes with membrane potential in one of the graphs
conductance_K(i) = g_K*(n(i)^4);
conductance_Na(i)=g_Na*(m(i)^3)*h(i);
% retrieving ionic currents
I_Na(i) = conductance_Na(i)*(V(i)-E_Na);
I_K(i) = conductance_K(i)*(V(i)-E_K);
I_Leak(i) = g_Leak*(V(i)-E_Leak);
% Calculating the input
Input = I_current(i) - (I_Na(i) + I_K(i) + I_Leak(i));
% Calculating the new membrane potential
V(i+1) = V(i) + Input* dt*(1/C);
% getting new values for the gating variables
m(i+1) = m(i) + (alphaM *(1-m(i)) - betaM * m(i))*dt;
n(i+1) = n(i) + (alphaN *(1-n(i)) - betaN * n(i))*dt;
h(i+1) = h(i) + (alphaH *(1-h(i)) - betaH * h(i))*dt;
end

figure('Name', 'Membrane Potential vs input')
plot(t(1/dt:end),V(0.01/dt:end-1), 'LineWidth', 2)
xlabel('Time (ms)')
end

% calculating alpha m and beta m based on known equations
function [alpha_m, beta_m] = m_equations(V, Vrest)
alpha_m = (2.5-0.1*(V-Vrest))/(exp(2.5-0.1*(V-Vrest))-1);
beta_m = 4*exp((Vrest-V)/18);
end
% calculate alpha n and beta n based on known equations
function [alpha_n, beta_n] = n_equations(V, Vrest)
alpha_n = (0.1-0.01*(V-Vrest))/(exp(1-0.1*(V-Vrest))-1);
beta_n = 0.125*exp((Vrest-V)/80);
end
% calculate alpha h and beta h based on known equations
function [alpha_h, beta_h] = h_equations(V, Vrest)
alpha_h = 0.07*exp((Vrest-V)/20);
beta_h = 1/(1+exp(3-0.1*(V-Vrest)));
end




