clc
clear all
close all
% Eric Walker 
% M.S. thesis NLLS

%% He, et al

load(['D:\Documents and Settings\Eric\My Documents\spring 2013\NASA'...
    ' Ames Data\B0005.mat']);
load(['D:\Documents and Settings\Eric\My Documents\spring 2013\NASA'...
    ' Ames Data\B0007.mat']);



first_batt = (-9.86e-7) * exp(5.752e-2 * (1:200)) + (8.983e-1) * ...
    exp((-8.340e-4) * (1:200)) + 0.005*randn(1,200);

second_batt = (-9.86e-7) * exp(5.752e-2 * (1:250)) + (8.983e-1) *...
    exp((-8.340e-4) * (1:250)) + 0.005*randn(1,250);



% Present the data in Figure 1.
hold on % Show all plots on the same figure.
plot(1:length(first_batt), first_batt, 'ko')
plot([1,200],first_batt(1)*0.8*[1,1],'k-','linewidth',1.5)
text(25,first_batt(1)*0.81,'EUL failure threshold')
xlabel('Cycle number (k)')
ylabel('Capacity (Ah)')
ylim([0.65, 0.91])
axis square 
box on

%NLLS
theta=[-9.86e-7,5.752e-2,8.983e-1,-8.34e-4];   


[theta_50_one_st, resnorm] = lsqnonlin(@(t_1) (second_batt(1:50))...
    - (-9.86e-7*exp(5.752e-2*...
    (1:50)) + t_1 * exp(-8.34e-4*(1:50))),8.983e-1); 

[theta_50_two_st, resnorm] = lsqnonlin(@(theta) (second_batt(1:50))...
    - (-9.86e-7*exp(5.752e-2*...
    (1:50)) + theta(1) * exp(theta(2)*(1:50))),[8.983e-1,-8.34e-4],...
    [0, -inf], [inf, 0]);

[theta_50_four_st, resnorm] = lsqnonlin(@(t_4) (second_batt(1:50)) ...
    - (t_4(1)*exp(t_4(2)*...
    (1:50)) + t_4(3) * exp(t_4(4)*(1:50))),[-9.86e-7,5.752e-2,...
    8.983e-1,-8.34e-4],[-inf, 0, 0, -inf], [0, inf, inf, 0]); 

[theta_50_four_st_const, resnorm] = lsqnonlin(@(t_4_c) ...
    (second_batt(1:50)) - (t_4_c(1)*exp(t_4_c(2)*...
    (1:50)) + t_4_c(3) * exp(t_4_c(4)*(1:50))),[-9.86e-7,...
    5.752e-2,8.983e-1,-8.34e-4],[-9.86e-7*1.05, 5.752e-2*0.95,...
    8.983e-1*0.95, -8.34e-4*1.05], [-9.86e-7*0.95, 5.752e-2*1.05,...
    8.983e-1*1.05, -8.34e-4*0.95]); 

[states,n_iter] = my_lsqnonlin_He(second_batt);

figure  % See three NLLS predictions
hold on
plot(1:200, (theta(1)*exp(theta(2)*...
    (1:200)) + theta_50_one_st * exp(theta(4)*(1:200))),'k-',...
    'linewidth',2)
plot(1:length(second_batt), (theta(1)*exp(theta(2)*...
    (1:length(second_batt))) + theta_50_two_st(1) * exp(...
    theta_50_two_st(2)*(1:length(second_batt)))),'k--',...
    'linewidth',2)
plot(1:length(second_batt), (theta_50_four_st(1)*exp(...
    theta_50_four_st(2)*...
    (1:length(second_batt))) + theta_50_four_st(3) * exp(...
    theta_50_four_st(4)*(1:length(second_batt)))),'k:',...
    'linewidth',2)
plot(1:length(second_batt), (theta_50_four_st_const(1)*exp(...
    theta_50_four_st_const(2)*...
    (1:length(second_batt))) + theta_50_four_st_const(3) * exp(...
    theta_50_four_st_const(4)*(1:length(second_batt)))),'k-.',...
    'linewidth',2)
plot(1:length(second_batt), (states(1)*exp(...
    states(2)*...
    (1:length(second_batt))) + states(3) * exp(...
    states(4)*(1:length(second_batt)))),'k.',...
    'linewidth',2)
plot(1:length(second_batt), second_batt, 'ko')
plot([1,250],second_batt(1)*0.8*[1,1],'k-','linewidth',1.5) 
text(25,second_batt(1)*0.81,'EUL failure threshold')
xlabel('Cycle number (k)')
ylabel('Capacity (Ah)')
ylim([0.65, 0.91])
axis square 
box on
legend('c state tracked ','c and d state tracked',...
    'four states tracked','four states tracked 5% constrained',...
    'NLLS code not lsqnonlin','observations') 

%% ECM

load(['D:\Documents and Settings\Eric\My Documents\spring 2013\'...
    'NASA Ames Data\B0006.mat']);

global first_discharge_time first_discharge_current...
    first_discharge_voltage;
first_discharge_voltage = B0006.cycle(1,2).data.Voltage_measured...
    (3:end);
first_discharge_time    = B0006.cycle(1,2).data.Time(3:end);
first_discharge_current = B0006.cycle(1,2).data.Current_measured...
    (3:end);


[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, [1.17e-8, 2.1,...
    1795.6, 0.28],[0 0 0 0], [1, 10, 2000, 1]);  
second_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured...
    (3:end);
second_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
second_discharge_current = B0006.cycle(1,4).data.Current_measured...
    (3:end);
second_discharge_time    = [second_discharge_time ...
    (second_discharge_time(1:100) + second_discharge_time(end))];
second_discharge_current = [second_discharge_current ...
    second_discharge_current(1:100)];


% Take the first forty data points for the NLLS prediction
first_discharge_voltage  = second_discharge_voltage(1:50);
first_discharge_time     = second_discharge_time(1:50);
first_discharge_current  = second_discharge_current(1:50);
[ecm_2, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, ecm,...
    [0 0 0 0], [1, 10, 2000, 1]); 


R = ecm_2(1);
Q = ecm_2(2);
C = ecm_2(3);
R_ct = ecm_2(4);

SOC_cell = 1 + second_discharge_current ./ (Q*3600) .*...
    second_discharge_time;

SOC_n = 0.79.*SOC_cell + 0.01;
SOC_p = 0.97-0.51*SOC_cell;

x_nsurf = SOC_n;
x_psurf_set = SOC_p;

U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 -...
    .0172./x_nsurf + ...
    .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984...
    * exp (...
    0.4465*x_nsurf - 0.4108);

U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 *...
    x_psurf_set.^4 + 342.909 * ...
    x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 *...
    x_psurf_set.^10 ) ./ ...
    ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 +...
    37.311 * x_psurf_set.^6 ...
    - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );

V_o = U_p_set - U_n ;

V_cell = V_o + second_discharge_current*R + Q./C .* exp(-...
    second_discharge_time./(R_ct * C))...
    + second_discharge_current.*R_ct.*(1-exp(-...
    second_discharge_time./(R_ct * C)));

% Now the constrained

[ecm_2_c, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, ecm,...
    0.95*ecm, 1.05*ecm); 


R = ecm_2_c(1);
Q = ecm_2_c(2);
C = ecm_2_c(3);
R_ct = ecm_2_c(4);

SOC_cell = 1 + second_discharge_current ./ (Q*3600) .*...
    second_discharge_time;

SOC_n = 0.79.*SOC_cell + 0.01;
SOC_p = 0.97-0.51*SOC_cell;

x_nsurf = SOC_n;
x_psurf_set = SOC_p;

U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 -...
    .0172./x_nsurf + ...
    .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984...
    * exp (...
    0.4465*x_nsurf - 0.4108);

U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 *...
    x_psurf_set.^4 + 342.909 * ...
    x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 *...
    x_psurf_set.^10 ) ./ ...
    ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 +...
    37.311 * x_psurf_set.^6 ...
    - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );

V_o = U_p_set - U_n ;

V_cell_c = V_o + second_discharge_current*R + Q./C .* exp(-...
    second_discharge_time./(R_ct * C))...
    + second_discharge_current.*R_ct.*(1-exp(-...
    second_discharge_time./(R_ct * C)));

first_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured...
    (3:end);
first_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
for y = 1:length(V_cell)
    if V_cell(y) < 2.5
        break
    end
end
for z = 1:length(V_cell_c)
    if V_cell_c(z) < 2.5
        break
    end
end
V_cell = V_cell(1:y);
V_cell_c = V_cell_c(1:z);
figure
hold on
plot(second_discharge_time(1:length(V_cell)),V_cell,'k-',...
    'linewidth',2)
plot(second_discharge_time(1:length(V_cell_c)),V_cell_c,'k--',...
    'linewidth',2)
plot(first_discharge_time,first_discharge_voltage, 'ko',...
    'linewidth', 1.5)
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis([0 4250 2.4 4])
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('NLLS prediction 928.6 (s)', ...
    'NLLS prediction 928.6 (s) 5% constrained', 'observations')

% Reset the following three variables after the above nlls.
first_discharge_voltage = B0006.cycle(1,2).data.Voltage_measured(3:end);
first_discharge_time    = B0006.cycle(1,2).data.Time(3:end);
first_discharge_current = B0006.cycle(1,2).data.Current_measured(3:end);
[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, ...
    [1.17e-8, 2.1, 1795.6, 0.28],[0 0 0 0], [1, 2.5, 2000, 1]);

%%  NLLS Single Particle model

[S_x_avg,resnorm,res] = lsqnonlin(@SP_obj_fun, ...
    [0.2447, 0.2536, 0.8957971, 0.5075848]...
    ,[],[],[],second_discharge_time(1:50),...
    second_discharge_current(1:50),second_discharge_voltage(1:50));

[S_x_avg_c,resnorm_c,res_c] = lsqnonlin(@SP_obj_fun, ...
    [0.2447, 0.2536, 0.8957971, 0.5075848]...
    ,[0.2447, 0.2536, 0.8957971, 0.5075848]*0.95,...
    [0.2447, 0.2536, 0.8957971, 0.5075848]*1.05,...
    [],second_discharge_time(1:50),...
    second_discharge_current(1:50),second_discharge_voltage(1:50));

[voltagePredi]        = SP(S_x_avg,second_discharge_time,...
    second_discharge_current);

[voltagePredi_c]      = SP(S_x_avg_c,second_discharge_time,...
    second_discharge_current);

for z = 1:length(voltagePredi_c)
    if voltagePredi_c(z) < 2.5
        break
    end
end
voltagePredi_c = voltagePredi_c(1:z);

figure
axis([0 4250 2.4 4])
axis square
hold on
xlabel('Time(s)');
ylabel('Voltage(V)');
plot(second_discharge_time(1:98),voltagePredi(1:98), 'k-',...
    'linewidth',2)
plot(second_discharge_time(1:length(voltagePredi_c)),voltagePredi_c, 'k--',...
    'linewidth',2)
plot(second_discharge_time(1:length(second_discharge_voltage))...
    ,second_discharge_voltage,'ko')
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
legend('NLLS prediction 928.6 (s)', ...
    'NLLS prediction 928.6 (s) 5% constrained', 'observations')