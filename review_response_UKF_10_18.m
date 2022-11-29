clc
clear all 
close all

%Eric Walker
%M.S. thesis UKF

%% UKF He, et al
second_batt = (-9.86e-7) * exp(5.752e-2 * (1:200)) + (8.983e-1) *...
    exp((-8.340e-4) * (1:200)) + 0.001*randn(1,200);

theta=[-9.86e-7; 5.752e-2; 8.983e-1; -8.34e-4];

P       = diag([(3.442e-8 - (-2.007e-6)), (6.221e-2 - 5.283e-2), ...
    (9.035e-1 - 8.931e-1), (-7.670e-4 - (-9.007e-4))])^2  ;...
    % 0.95 confidence from He, et al.
kappa   = 0.5;
Q_storage = []
tic
for i = 1:length(second_batt)
    C        = chol(P);
    Chi(:,1) = theta;
    Chi(:,2) = theta + sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,3) = theta + sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,4) = theta + sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,5) = theta + sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    Chi(:,6) = theta - sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,7) = theta - sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,8) = theta - sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,9) = theta - sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    W(1)     = kappa / sqrt(4 + kappa);
    W(2:9)   = 1     / (2*(4+kappa));
    W        = W/sum(W);

    Q_cap    = Chi(1,:) .* exp( Chi(2,:) * i ) + Chi(3,:) .* ...
        exp( Chi(4,:) * i ); 
    
    Q_hat    = sum(W.*Q_cap);
    Q_storage = [Q_storage; Q_hat];
    
    if i < 51
        P_yy     = sum(W.*((Q_cap).^2));
        
        
        P_xy     = (Chi - repmat(theta,1,9)) .* [W;W;W;W] *(Q_cap...
            - Q_hat)'; 
        K_k      = P_xy/(P_yy+0.05);
        egg      = theta + K_k*(second_batt(i) - Q_hat);
        theta(1) = egg(1);
        theta(2) = egg(2);
        theta(3) = egg(3);
        theta(4) = egg(4);


        P        = P - K_k * P_xy';
    end


end
disp('He et al UKF time')
toc

figure
hold on
plot(1:length(second_batt), Q_storage, 'k-','linewidth',2)
plot(1:length(second_batt), second_batt,'ko','linewidth',1.5)
legend('UKF prediction k=50', 'observations')
plot([1,200],second_batt(1)*0.8*[1,1],'k-','linewidth',1.5) 
text(25,second_batt(1)*0.81,'EUL failure threshold')
ylim([0.65, 0.91])
xlabel('k, Cycle index (cycle)')
ylabel('Q, Capacity (Ah)')
axis square
box on


%% UKF ECM
clear all
load(['B0006.mat']);


global first_discharge_time first_discharge_current...
    first_discharge_voltage;
first_discharge_voltage = B0006.cycle(1,2).data.Voltage_measured...
    (3:end);
first_discharge_time    = B0006.cycle(1,2).data.Time(3:end);
first_discharge_current = B0006.cycle(1,2).data.Current_measured...
    (3:end);

second_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured...
    (3:end);
second_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
second_discharge_current = B0006.cycle(1,4).data.Current_measured...
    (3:end);
second_discharge_time    = [second_discharge_time ...
    (second_discharge_time(1:30) + second_discharge_time(end))];
second_discharge_current = [second_discharge_current ...
    second_discharge_current(1:30)];

[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, ...
    [1.17e-8, 2.1, 1795.6, 0.28],[0 0 0 0], [1, 10, 2000, 2]);
ecm=ecm'
sigma = 0.0015;


V_cell_storage = [];
R = ecm(1) %+ ecm(1)/4 * randn(1,50)%1795;
Q = ecm(2) % + ecm(2)/20*(randn(1,50))); % %2.0593;
C = ecm(3) %+ ecm(3)/4 * randn(1,50);%1.17e-8;
R_ct = ecm(4) %+ ecm(4)/4*randn(1,50) %0.1451;


P       = (diag(ecm) / 20)^2;


ecm_orig= ecm;
kappa = 0.5;
tic
for i = 1:length(second_discharge_voltage)

    C        = chol(P);
    Chi(:,1) = ecm;
    Chi(:,2) = ecm + sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,3) = ecm + sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,4) = ecm + sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,5) = ecm + sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    Chi(:,6) = ecm - sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,7) = ecm - sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,8) = ecm - sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,9) = ecm - sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    W(1)     = kappa / sqrt(4 + kappa);
    W(2:9)   = 1     / (2*(4+kappa));
    W        = W/sum(W);

    SOC_cell = 1 + second_discharge_current(i) ./ ...
        (Chi(2,:)*3600) .* second_discharge_time(i);

    SOC_n = 0.79.*SOC_cell + 0.01;
    SOC_p = 0.97-0.51*SOC_cell;
    
    x_nsurf = SOC_n;
    x_psurf_set = SOC_p;

    
    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 -...
        .0172./x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -...
        .7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 *...
        x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 *...
        x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 +...
        37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    
    V_o = U_p_set - U_n ;
    
    V_cell = V_o + second_discharge_current(i)*Chi(1,:) + ...
        Chi(2,:)./Chi(3,:) .* exp(-second_discharge_time(i)./...
        (Chi(4,:) .* Chi(3,:)))...
        + second_discharge_current(i).*Chi(4,:).*(1-exp(...
        -second_discharge_time(i)./(Chi(4,:) .* Chi(3,:))));

    
    y_hat    = sum(W.*V_cell);
    V_cell_storage = [V_cell_storage; y_hat];
    
    if i < 50
        P_yy     = sum(W.*((V_cell).^2));
        
        
        P_xy     = (Chi - repmat(ecm,1,9)) .* [W; W; W; W] *...
            (V_cell - y_hat)'; %P_xz is (1x1).
        K_k      = P_xy/(P_yy+0.05);
        ecm      = ecm + K_k*(second_discharge_voltage(i) - y_hat);
        
        
        err(i)   = (second_discharge_voltage(i) - y_hat);

        P        = P - K_k * P_xy';
    end


end

% This block of code is to stop the prediction when the
%voltage drops below 2.5 volts.
k=1
while V_cell_storage(k) > 2.5 && k < 195
    k = k+1;
end
V_cell_storage = V_cell_storage(1:k);
disp('ECM UKF time')
toc

figure
hold on
plot(second_discharge_time(1:length(V_cell_storage)),...
    V_cell_storage,'k-','linewidth',2)
plot(second_discharge_time(1:length(second_discharge_voltage)),...
    second_discharge_voltage, 'ko', 'linewidth', 1.5)
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis([0 4250 2.4 4])
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('UKF prediction 928.6 (s)', 'observations')

%% UKF SP

S_n    = 0.2604;
S_p    = 0.2570;
k_n    = 37.4312e-12;
k_p    = 17.4733e-12;
R_n    = 2e-6;
R_p    = 2e-6;
D_n    = 29.0798e-15;
D_p    = 27.9034e-15;
c_nmax = 30074.5;
c_pmax = 51563.5;
c_e    = 1000;
x_navg = 0.9401;
x_pavg = 0.5169;
T      = 298.15;
R_g    = 8.3143;
F      = 96487;
alpha_a  = 0.5;
alpha_c  = 0.5;

P = (diag([S_n,S_p,x_navg,x_pavg]./20))^2;
C = P;
kappa = 0.2;
tic
for i=1:length(second_discharge_time);


    C        = chol(P);
    Chi(:,1) = [S_n; S_p; x_navg; x_pavg];
    Chi(:,2) = [S_n; S_p; x_navg; x_pavg] + ...
        sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,3) = [S_n; S_p; x_navg; x_pavg] + ...
        sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,4) = [S_n; S_p; x_navg; x_pavg] + ...
        sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,5) = [S_n; S_p; x_navg; x_pavg] + ...
        sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    Chi(:,6) = [S_n; S_p; x_navg; x_pavg] - ...
        sqrt(4 + kappa)*[C(1,1); 0; 0; 0];
    Chi(:,7) = [S_n; S_p; x_navg; x_pavg] - ...
        sqrt(4 + kappa)*[0; C(2,2); 0; 0];  
    Chi(:,8) = [S_n; S_p; x_navg; x_pavg] - ...
        sqrt(4 + kappa)*[0; 0; C(3,3); 0];  
    Chi(:,9) = [S_n; S_p; x_navg; x_pavg] - ...
        sqrt(4 + kappa)*[0; 0; 0; C(4,4)];  
    W(1)     = kappa / sqrt(4 + kappa);
    W(2:9)   = 1     / (2*(4+kappa));
    W        = W/sum(W); 



    %%%%%Now the SP measurement model
    Iapp = second_discharge_current(i);
    J_n  = -Iapp./Chi(1,:);
    J_p  = Iapp./Chi(2,:);
    x_nsurf = Chi(3,:) - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
    x_psurf = Chi(4,:) - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);

    U_n     = .8214 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./...
        x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 ...
        * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p     = ( -4.8801 + 88.669 * x_psurf.^2 - 401.119 * x_psurf.^4 ...
    + 342.909 * ...
        x_psurf.^6 - 462.471 * x_psurf.^8 + 433.434 * x_psurf.^10 )...
        ./ ...
        ( -1 + 18.933*x_psurf.^2 - 79.532 * x_psurf.^4 + 37.311 *...
        x_psurf.^6 ...
        - 73.083 * x_psurf.^8 + 95.96*x_psurf.^10 );
    eta_n   = R_g * T ./ (F * alpha_a) .* log( ( J_n + (-4*c_e*...
        F.^2*c_nmax.^2*k_n.^2.*x_nsurf.^2 ...
       + 4*c_e*F^2*c_nmax.^2*k_n.^2.*x_nsurf+J_n.^2).^0.5 ) ./ ...
       (2*F*c_e^0.5*k_n.*(c_nmax.*x_nsurf).^0.5 .* ...
       (c_nmax-c_nmax.*x_nsurf).^0.5) );
   
    eta_p   = R_g * T / (F * alpha_c) .* log( ( J_p + (-4*c_e...
        *F^2*c_pmax^2*k_p^2.*x_psurf.^2 ...
       + 4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf+J_p.^2).^0.5 ) ./ ...
       (2*F*c_e^0.5*k_p.*(c_pmax.*x_psurf).^0.5 .* ...
       (c_pmax-c_pmax.*x_psurf).^0.5) );
   
    zeta = U_p + eta_p - U_n - eta_n;
    y_hat    = sum(W.*zeta);
    P_yy     = sum(W.*((zeta-y_hat).^2));
    y_hat_storage(i) = y_hat;
    
    if i <= 50
        P_xy     = (Chi - repmat([S_n; S_p; x_navg; x_pavg],1,9))...
            .* [W; W; W; W] *(zeta - y_hat)'; 
        K_k      = P_xy/(P_yy+0.2);
        egg      = [S_n; S_p; x_navg; x_pavg] + K_k*(...
            second_discharge_voltage(i) - y_hat);
        S_n      = egg(1);
        S_p      = egg(2);
        x_navg   = egg(3);
        x_pavg   = egg(4);
        
       Q        = 0.4e-4 * ones(4,4);
       P        = P - K_k * P_xy' + Q;
    end
        J_n  = -Iapp./S_n;
        J_p  = Iapp./S_p;
        x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
        x_pavg    = x_pavg - 3 * J_p / (F * R_p * c_pmax) ;
        
    
end
disp('SP UKF computation time')
toc

y_hat_storage = y_hat_storage(1:find(y_hat_storage<2.5,1));
figure
hold on
plot(second_discharge_time(1:length(y_hat_storage)), y_hat_storage,...
    'k-', 'linewidth', 2)
axis([0 4300 2.4 4])
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('UKF prediction 928.6 (s)', 'observations')
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
plot(second_discharge_time(1:length(second_discharge_voltage)),...
    second_discharge_voltage,'ko')
axis square
box on