clc
clear all
close all
% Eric Walker 
% M.S. thesis

%% 2. He, et al model

%Load the data set.

load('D:\Documents and Settings\Eric\My Documents\spring 2013\NASA Ames Data\B0005.mat');
load('D:\Documents and Settings\Eric\My Documents\spring 2013\NASA Ames Data\B0007.mat');



first_batt = (-9.86e-7) * exp(5.752e-2 * (1:200)) + (8.983e-1) * exp((-8.340e-4) * (1:200)) + 0.005*randn(1,200);

second_batt = (-9.86e-7) * exp(5.752e-2 * (1:200)) + (8.983e-1) * exp((-8.340e-4) * (1:200)) + 0.005*randn(1,200);



% Present the data in Figure 1.
hold on % Show all plots on the same figure.
plot(1:length(first_batt), first_batt, 'ko')
plot([1,200],first_batt(1)*0.8*[1,1],'k-','linewidth',1.5) % The EUL line.
text(25,first_batt(1)*0.81,'EUL failure threshold')
xlabel('Cycle number (k)')
ylabel('Capacity (Ah)')
ylim([0.65, 0.91])
axis square 
box on

%NLLS
theta=[-9.86e-7,5.752e-2,8.983e-1,-8.34e-4];  % A similar battery is tested 


[theta_50, resnorm] = lsqnonlin(@(theta) (second_batt(1:50)) - (theta(1)*exp(theta(2)*...
    (1:50)) + theta(3) * exp(theta(4)*(1:50))),zeros(1,4)); 

[theta_100, resnorm] = lsqnonlin(@(theta) (second_batt(1:100)) - (theta(1)*exp(theta(2)*...
    (1:100)) + theta(3) * exp(theta(4)*(1:100))),zeros(1,4)); 

[theta_150, resnorm] = lsqnonlin(@(theta) (second_batt(1:150)) - (theta(1)*exp(theta(2)*...
    (1:150)) + theta(3) * exp(theta(4)*(1:150))),zeros(1,4));

figure  % See three NLLS predictions
hold on
plot(1:300, (theta_50(1)*exp(theta_50(2)*...
    (1:300)) + theta_50(3) * exp(theta_50(4)*(1:300))),'k-','linewidth',2)
plot(1:length(second_batt), (theta_100(1)*exp(theta_100(2)*...
    (1:length(second_batt))) + theta_100(3) * exp(theta_100(4)*(1:length(second_batt)))),'k--','linewidth',2)
plot(1:length(second_batt), (theta_150(1)*exp(theta_150(2)*...
    (1:length(second_batt))) + theta_150(3) * exp(theta_150(4)*(1:length(second_batt)))),'k-.','linewidth',2)
plot(1:length(second_batt), second_batt, 'ko')
plot([1,200],second_batt(1)*0.8*[1,1],'k-','linewidth',1.5) % The EUL line.
text(25,second_batt(1)*0.81,'EUL failure threshold')
xlabel('Cycle number (k)')
ylabel('Capacity (Ah)')
ylim([0.65, 0.91])
axis square 
box on

legend('RUL prediction k=50','RUL prediction k=100','RUL prediction k=150','observations') 

%UKF
theta=[-9.86e-7; 5.752e-2; 8.983e-1; -8.34e-4];

P       = diag([(3.442e-8 - (-2.007e-6)), (6.221e-2 - 5.283e-2), ...
    (9.035e-1 - 8.931e-1), (-7.670e-4 - (-9.007e-4))])^2  ; % 0.95 confidence from He, et al.
kappa   = 0.5;
Q_storage = []
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
    W(1)     = kappa / sqrt(1 + kappa);
    W(2:9)   = 1     / (2*(1+kappa));
    W        = W/sum(W);

    Q_cap    = Chi(1,:) .* exp( Chi(2,:) * i ) + Chi(3,:) .* exp( Chi(4,:) * i ); 
    
    Q_hat    = sum(W.*Q_cap);
    Q_storage = [Q_storage; Q_hat];
    
    if i < 50
        P_yy     = sum(W.*((Q_cap).^2));
        
        
        P_xy     = (Chi - repmat(theta,1,9)) .* [W;W;W;W] *(Q_cap - Q_hat)'; %P_xz is (1x1).
        K_k      = P_xy/(P_yy+0.05);
        egg      = theta + K_k*(second_batt(i) - Q_hat);
        theta(1) = egg(1);
        theta(2) = egg(2);
        theta(3) = egg(3);
        theta(4) = egg(4);


        P        = P - K_k * P_xy';
    end


end
figure
hold on
plot(1:length(second_batt), Q_storage, 'k-','linewidth',2)
plot(1:length(second_batt), second_batt,'ko','linewidth',1.5)
legend('RUL prognosis k=50', 'observations')
plot([1,200],second_batt(1)*0.8*[1,1],'k-','linewidth',1.5) % The EUL line.
text(25,second_batt(1)*0.81,'EUL failure threshold')
ylim([0.65, 0.91])
xlabel('Cycle index (k)')
ylabel('Capacity (Ah)')
axis square
box on



%Now particle filter B0007.  a and c, together, are better than b and d,
%together.
theta_set=repmat(theta,1,200);
theta_set(1,1:50) = theta(1)% + 3*theta(1) * (0.5-rand(50,1));
theta_set(2,51:100) = theta(2)% + theta(2)/2 * randn(50,1);
theta_set(3,101:150) = theta(3) + theta(3)/10 * (0.5-rand(50,1));
theta_set(4,151:200) = theta(4)% + theta(4)/3 * randn(50,1);
weights = 0.005 * ones(1,200);

for j = 1:200
    choose_par(j,:) = theta_set(1,j) * exp(theta_set(2,j) * (1:200)) + theta_set(3,j) * exp(theta_set(4,j)*(1:200));
end






% Make the RUL
for i = 1:200
RULs(i) = fsolve(@(k) 0.8*(second_batt(1)) - (theta_set(1,i)*exp(theta_set(2,i)*...
    k) + theta_set(3,i) * exp(theta_set(4,i)*k)), 180);
end

sigma = 0.1;
for i = 1:200
    if i ==50
        weights_50 = weights;
    end
    if i == 100
        weights_100 = weights;
    end
    if i == 150
        weights_150 = weights;
    end
    % Get the likelihood
    likelihood = 1/(sigma*sqrt(2*pi)) * exp(-1/2 * ((second_batt(i)) - ...
        (theta_set(1,:) .* exp(theta_set(2,:) * i) + theta_set(3,:) .* exp(theta_set(4,:) * i))).^2 / sigma^2);
    % Update the weights
    weights = weights .* likelihood;
    weights = weights / sum(weights);
end

[RULs, ind] = sort(RULs);
weights_50s = weights_50(ind);
weights_100s = weights_100(ind);
weights_150s = weights_150(ind);

figure

xlabel('Cycle index (k)')
ylabel('Capacity (Ah)')
axis square
hold on 
plot(RULs', weights_50s*10 + 0.8*second_batt(1),'k-', 'linewidth', 2)
plot(RULs', weights_100s*10 + 0.8*second_batt(1),'k--', 'linewidth', 2)
plot(RULs', weights_150s*10 + 0.8*second_batt(1),'k:', 'linewidth', 2)
plot(1:length(second_batt), second_batt,'ko','linewidth',1.5)


% Make them range from 0 to 1, otherwise they will be light.
for j = 1:200
    plot(1:50,choose_par(j,1:50),'color',(1-weights_50(j)/max(weights_50))*[1, 1, 1]); % Smaller numbers are darker.
    plot(50:100,choose_par(j,50:100),'color',(1-weights_100(j)/max(weights_100))*[1, 1, 1]); % Smaller numbers are darker.
    plot(100:150,choose_par(j,100:150),'color',(1-weights_150(j)/max(weights_150))*[1, 1, 1]); % Smaller numbers are darker.
    plot(150:200,choose_par(j,150:200), 'color',(1-weights(j)/max(weights))*[1, 1, 1]);
end

plot(RULs', weights_50s*10 + 0.8*second_batt(1),'k-', 'linewidth', 2)
plot(RULs', weights_100s*10 + 0.8*second_batt(1),'k--', 'linewidth', 2)
plot(RULs', weights_150s*10 + 0.8*second_batt(1),'k:', 'linewidth', 2)
plot(1:length(second_batt), second_batt,'ko','linewidth',1.5)
legend('RUL prognosis k=50', 'k=100', 'k=150', 'observations')
plot([1,200],second_batt(1)*0.8*[1,1],'k-','linewidth',1.5) % The EUL line.
text(25,second_batt(1)*0.81,'EUL failure threshold')
ylim([0.65, 0.91])
xlabel('Cycle index (k)')
ylabel('Capacity (Ah)')
axis square
box on

err_early = sum(weights_50s.*RULs)-190
err_late  = sum(weights_100s.*RULs)-190
err_final = sum(weights_150s.*RULs)-190
sig_early = sqrt(sum(weights_50s.*(RULs - (err_early + 190)).^2) )
sig_late  = sqrt(sum(weights_100s.*(RULs - (err_late + 190)).^2) )
sig_final = sqrt(sum(weights_150s.*(RULs - (err_final + 190)).^2) )


%% 3. Equivalent circuit model
clear all
load('D:\Documents and Settings\Eric\My Documents\spring 2013\NASA Ames Data\B0006.mat');

global first_discharge_time first_discharge_current first_discharge_voltage;
first_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured(3:end);
first_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
first_discharge_current = B0006.cycle(1,4).data.Current_measured(3:end);


[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, [1.17e-8, 2.1, 1795.6, 0.28],[0 0 0 0], [2000, 2.5, 1, 1]);  % Maybe include lb and ub , [0 0 0 0], [2000, 2.5, 1, 1] .
second_discharge_voltage = B0006.cycle(1,2).data.Voltage_measured(3:end);
second_discharge_time    = B0006.cycle(1,2).data.Time(3:end);
second_discharge_current = B0006.cycle(1,2).data.Current_measured(3:end);
second_discharge_time    = [second_discharge_time (second_discharge_time(1:30) + second_discharge_time(end))];
second_discharge_current = [second_discharge_current second_discharge_current(1:30)];

figure
hold on
plot(first_discharge_time,first_discharge_voltage, 'ko', 'linewidth', 1.5)
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis([0 4250 2.4 4])
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')

% Take the first forty data points for the NLLS prediction
first_discharge_voltage  = second_discharge_voltage(1:50);
first_discharge_time     = second_discharge_time(1:50);
first_discharge_current  = second_discharge_current(1:50);
[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, [1.17e-8, 2.1, 1795.6, 0.28],[0 0 0 0], [2000, 2.5, 1, 1]); % R Q C R_ct



R = ecm(1);
Q = ecm(2);
C = ecm(3);
R_ct = ecm(4);

SOC_cell = 1 + second_discharge_current ./ (Q*3600) .* second_discharge_time;

SOC_n = 0.79.*SOC_cell + 0.01;
SOC_p = 0.97-0.51*SOC_cell;

x_nsurf = SOC_n;
x_psurf_set = SOC_p;

U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
    .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
    0.4465*x_nsurf - 0.4108);

U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
    x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
    ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
    - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );

V_o = U_p_set - U_n ;

V_cell = V_o + second_discharge_current*R + Q./C .* exp(-second_discharge_time./(R_ct * C))...
    + second_discharge_current.*R_ct.*(1-exp(-second_discharge_time./(R_ct * C)));

first_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured(3:end);
first_discharge_time    = B0006.cycle(1,4).data.Time(3:end);

figure
hold on
plot(first_discharge_time,first_discharge_voltage, 'ko', 'linewidth', 1.5)
plot(second_discharge_time,V_cell,'k-','linewidth',2)
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis([0 4250 2.4 4])
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Observations (V)', 'NLLS prediction 928.6 (s)')

% Reset the following three variables after the above nlls.
first_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured(3:end);
first_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
first_discharge_current = B0006.cycle(1,4).data.Current_measured(3:end);
[ecm, resnorm, residuals] = lsqnonlin( @ecm_obj_fun, [1.17e-8, 2.1, 1795.6, 0.28],[0 0 0 0], [2000, 2.5, 1, 1]);

second_discharge_time    = [second_discharge_time (second_discharge_time(1:30) + second_discharge_time(end))];
second_discharge_current = [second_discharge_current second_discharge_current(1:30)];

% UKF




ecm
% residuals
sigma = 0.0015;
% param_set = repmat([ecm(1); ecm(2)],1,100);
% param_set(1,1:50) = ecm(1) + ecm(1)/4*randn(50,1);
% param_set(2,51:100)  = ecm(2) + ecm(2)/4*randn(50,1);
% weights = 0.05*ones(1,100);


V_cell_storage = [];
R = ecm(1) %+ ecm(1)/4 * randn(1,50)%1795;
Q = (ecm(2)) % + ecm(2)/20*(randn(1,50))); % %2.0593;
C = ecm(3) %+ ecm(3)/4 * randn(1,50);%1.17e-8;
R_ct = ecm(4) %+ ecm(4)/4*randn(1,50) %0.1451;

weights = 0.02*ones(1,50);
P       = ecm(2) / 20;
kappa = 0.5;
for i = 1:length(second_discharge_voltage)
    C        = chol(P);
    Chi(:,1) = Q;
    Chi(:,2) = Q + sqrt(1 + kappa)*[C(1,1)];
    Chi(:,3) = Q - sqrt(1 + kappa)*[C(1,1)];  
     
    W(1)     = kappa / sqrt(1 + kappa);
    W(2:3)   = 1     / (2*(1+kappa));
    W        = W/sum(W);

    SOC_cell = 1 + second_discharge_current(i) ./ (Chi*3600) .* second_discharge_time(i);

    SOC_n = 0.79.*SOC_cell + 0.01;
    SOC_p = 0.97-0.51*SOC_cell;
    
    x_nsurf = SOC_n;
    x_psurf_set = SOC_p;

    
    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    
    V_o = U_p_set - U_n ;
    
    V_cell = V_o + second_discharge_current(i)*R + Q./C .* exp(-second_discharge_time(i)./(R_ct * C))...
        + second_discharge_current(i).*R_ct.*(1-exp(-second_discharge_time(i)./(R_ct * C)));

    
    y_hat    = sum(W.*V_cell);
    V_cell_storage = [V_cell_storage; y_hat];
    
    if i < 50
        P_yy     = sum(W.*((V_cell).^2));
        
        
        P_xy     = (Chi - repmat([Q],1,3)) .* [W] *(V_cell - y_hat)'; %P_xz is (1x1).
        K_k      = P_xy/(P_yy+0.05);
        Q      = [Q] + K_k*(second_discharge_voltage(i) - y_hat);
        
        
        err(i)   = (second_discharge_voltage(i) - y_hat);
        Q_noise        = 0.4e-4 ;
        P        = P - K_k * P_xy' + Q_noise;
    end


end

% This block of code is to stop the prediction when the voltage drops below
% 2.5 volts.
k=1
while V_cell_storage(k) > 2.5 && k < 195
    k = k+1;
end
V_cell_storage = V_cell_storage(1:k);
%
figure
hold on
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage, 'ko', 'linewidth', 1.5)
plot(second_discharge_time(1:length(V_cell_storage)),V_cell_storage,'k-','linewidth',2)
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis([0 4250 2.4 4])
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Observations (V)', 'UKF prediction 928.6 (s)')

% PF
sigma = 0.0015;
% param_set = repmat([ecm(1); ecm(2)],1,100);
% param_set(1,1:50) = ecm(1) + ecm(1)/4*randn(50,1);
% param_set(2,51:100)  = ecm(2) + ecm(2)/4*randn(50,1);
% weights = 0.05*ones(1,100);

V_cell_storage = [];
R = ecm(1) %+ ecm(1)/4 * randn(1,50)%1795;
Q = (ecm(2) + ecm(2)/12*(randn(1,50))); % %2.0593;
C = ecm(3) %+ ecm(3)/4 * randn(1,50);%1.17e-8;
R_ct = ecm(4) %+ ecm(4)/4*randn(1,50) %0.1451;


weights = 0.02*ones(1,50);
V_cell_storage = [];
for i = 1:length(second_discharge_time)
    tic
    % Get the likelihood
    % Quantity inside the square, first
    SOC_cell = 1 + second_discharge_current(i) ./ (Q*3600) .* second_discharge_time(i);

    SOC_n = 0.79.*SOC_cell + 0.01;
    SOC_p = 0.97-0.51*SOC_cell;
    
    x_nsurf = SOC_n;
    x_psurf_set = SOC_p;
    
    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    
    V_o = U_p_set - U_n ;
    
    V_cell = V_o + second_discharge_current(i)*R + Q./C .* exp(-second_discharge_time(i)./(R_ct * C))...
        + second_discharge_current(i).*R_ct.*(1-exp(-second_discharge_time(i)./(R_ct * C)));
    V_cell_storage = [V_cell_storage; V_cell];
  
    if i<=194
    quantity = (second_discharge_voltage(i) - V_cell).^2;
    likelihood = 1./(sigma*sqrt(2*pi)) .* exp(-1/2 * (quantity).^2 ./ sigma^2);
    % Update the weights
    weights = weights .* likelihood;
    weights = weights/sum(weights);
    end
    if i == 50
        weights_50 = weights;
    end
    if i == 100
        weights_100 = weights;
    end
    if i == 150
        weights_150 = weights;
    end
toc
end
disp('Now one-time time')



% figure
k=1
for j = 1:50
    vec = V_cell_storage(:,j);
    while vec(k,1)>2.5  && k<length(second_discharge_time)
        k=k+1;
    end
    EODs(j) = k;
    for_plot(j,1:k) = vec(1:k)';
    k=1;
end
toc

[EODs ind] = sort(EODs);
weights_50s = weights_50(ind);
weights_100s = weights_100(ind);
weights_150s = weights_150(ind);

figure
hold on
plot(second_discharge_time(EODs), weights_50s*5 + 2.5,'k-','linewidth',2)
plot(second_discharge_time(EODs), weights_100s*5 + 2.5,'k--','linewidth',2)
plot(second_discharge_time(EODs), weights_150s*5 + 2.5,'k:','linewidth',2)
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage, 'ko', 'linewidth', 1.5)

axis([0 4250 2.4 4])
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis square
box on
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('EOD prognosis 928.6 (s)','1,853.0 (s)', '2,802.3 (s)','observations')
for j=1:50
    plot(second_discharge_time(1:50), for_plot(j,1:50), 'color', (1-weights_50(j)/max(weights_50))*[1,1,1])
    plot(second_discharge_time(50:100), for_plot(j,50:100), 'color', (1-weights_100(j)/max(weights_100))*[1,1,1])
    plot(second_discharge_time(100:150), for_plot(j,100:150), 'color', (1-weights_150(j)/max(weights_150))*[1,1,1])
    plot(second_discharge_time(150:length(for_plot)), for_plot(j,150:end), 'color', (1-weights_150(j)/max(weights_150))*[1,1,1])
 
end
plot(second_discharge_time(EODs), weights_50s*5 + 2.5,'k-','linewidth',2)
plot(second_discharge_time(EODs), weights_100s*5 + 2.5,'k--','linewidth',2)
plot(second_discharge_time(EODs), weights_150s*5 + 2.5,'k:','linewidth',2)
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage, 'ko', 'linewidth', 1.5)
figure



err_early = sum(weights_50.*second_discharge_time(EODs))-3690
err_late  = sum(weights_100.*second_discharge_time(EODs))-3690
err_final = sum(weights_150.*second_discharge_time(EODs))-3690
sig_early = sqrt(sum(weights_50.*(second_discharge_time(EODs) - (err_early + 3690)).^2) )
sig_late  = sqrt(sum(weights_100.*(second_discharge_time(EODs) - (err_late + 3690)).^2) )
sig_final = sqrt(sum(weights_150.*(second_discharge_time(EODs) - (err_final + 3690)).^2) )
sum(weights_50.*second_discharge_time(EODs))
%%  4. Single Particle model
clc
clear all


%Load the data set.
load('D:\Documents and Settings\Eric\My Documents\spring 2013\NASA Ames Data\B0006.mat')

% When assigning the data set to variables, truncate the first two points, 
% because they are inconsistent with the rest of the data.

first_discharge_voltage = B0006.cycle(1,2).data.Voltage_measured(3:end);
first_discharge_time    = B0006.cycle(1,2).data.Time(3:end);
first_discharge_current = B0006.cycle(1,2).data.Current_measured(3:end);

   
[S_R_fit,resnorm,res] = lsqnonlin(@many_dof_obj_fun, [3.41,3.86,0], [0, 0, 0],[4, 4, 0.2],[],...
    first_discharge_time,first_discharge_current,first_discharge_voltage);
S_nfit  = S_R_fit(1)
S_pfit  = S_R_fit(2)
R       = S_R_fit(3)
% Now run the model over again, with S_nfit and S_pfit.

% V_cell_fit = first_discharge_voltage - res;
% 
% figure  % the fitted training plot
% plot(first_discharge_time,V_cell_fit,'k-', first_discharge_time, ...
%     first_discharge_voltage,'ko')
% xlabel('time (s)')
% ylabel('V_c_e_l_l (V)')

%%%%%%%%%%%%%%%%% NLLS section
% [S_R_fit,resnorm,res] = lsqnonlin(@many_dof_obj_fun, [3.41,3.86,0], [0, 0, 0],[4, 4, 0.2],[],...
%     first_discharge_time(1:50),first_discharge_current(1:50),first_discharge_voltage(1:50));
% S_n  = S_R_fit(1)
% S_p  = S_R_fit(2)
% R       = S_R_fit(3)
% 
% 
% k_n    = 37.4312e-12;
% k_p    = 17.4733e-12;
% R_n    = 2e-6;
% R_p    = 2e-6;
% D_n    = 29.0798e-15;
% D_p    = 27.9034e-15;
% c_nmax = 30074.5;
% c_pmax = 51563.5;
% c_e    = 1000;
% x_navg = 0.8957971;
% x_pavg_set = 0.5075848;
% T      = 298.15;
% R_g    = 8.3143;
% F      = 96487;
% alpha_a  = 0.5;
% alpha_c  = 0.5;
% 
% for i = 1:(length(first_discharge_voltage)+29)
% %     tic
%     Iapp = first_discharge_current(i);
%     J_n  = -Iapp./S_n;
%     J_p  = Iapp./S_p;
%     x_nsurf = x_navg - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
%     x_psurf_set = x_pavg_set - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);
% 
%     U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
%         .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
%         0.4465*x_nsurf - 0.4108);
%     
%     U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
%         x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
%         ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
%         - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
%     eta_n    = R_g * T ./ (F * alpha_a) .* log( ( J_n + (-4*c_e*F.^2*c_nmax.^2*k_n.^2.*x_nsurf.^2 ...
%        + 4*c_e*F^2*c_nmax.^2*k_n.^2.*x_nsurf+J_n.^2).^0.5 ) ./ (2*F*c_e^0.5*k_n.*(c_nmax.*x_nsurf).^0.5 .* ...
%        (c_nmax-c_nmax.*x_nsurf).^0.5) ) - Iapp./S_n .* R_film;
%    
%     eta_p_set    = R_g * T / (F * alpha_c) .* log( ( J_p + (-4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set.^2 ...
%        + 4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set+J_p.^2).^0.5 ) ./ (2*F*c_e^0.5*k_p.*(c_pmax.*x_psurf_set).^0.5 .* ...
%        (c_pmax-c_pmax.*x_psurf_set).^0.5) );
%    
%     V_cell_set(i,:) = real(U_p_set + eta_p_set - U_n - eta_n);
% 
%     x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
%     x_pavg_set    = x_pavg_set - 3 * J_p / (F * R_p * c_pmax) ;
% end
% figure
% hold on
% plot(first_discharge_time,first_discharge_voltage,'ko')
% plot(first_discharge_time,V_cell_set,'k-','linewidth',2)
% figure

%%%%%%%%%%%%%%%%%

% Now compare PF and UKF 
S_n    = 0.2447;
S_p    = 0.2536;
k_n    = 37.4312e-12;
k_p    = 17.4733e-12;
R_n    = 2e-6;
R_p    = 2e-6;
D_n    = 29.0798e-15;
D_p    = 27.9034e-15;
c_nmax = 30074.5;
c_pmax = 51563.5;
c_e    = 1000;
x_navg = 0.8957971;
x_pavg_set = 0.5075848;
T      = 298.15;
R_g    = 8.3143;
F      = 96487;
alpha_a  = 0.5;
alpha_c  = 0.5;

%S_p  = 0.2536 + 0.02*randn(50,1);  % PF is making uncertainty
% about the state of charge, x_pavg_set.
%S_n    = 0.2447 + 0.05*randn(50,1);
%S_p     =  [(0.2536 - 0.003*sort(randn(50,1)))];
R_film  = [( 0.01 * sort(randn(50,1)))]
x_navg = 0.8957971 - 0.035*sort(randn(50,1));
weights    = 0.02*ones(1,50);

second_discharge_voltage = B0006.cycle(1,4).data.Voltage_measured(3:end);
second_discharge_time    = B0006.cycle(1,4).data.Time(3:end);
second_discharge_current = B0006.cycle(1,4).data.Current_measured(3:end);
second_discharge_time    = [second_discharge_time (second_discharge_time(1:30) + second_discharge_time(end))];
second_discharge_current = [second_discharge_current second_discharge_current(1:30)];
stop_cycle = length(second_discharge_voltage)+29;
tic
for i = 1:stop_cycle
%     tic
    Iapp = second_discharge_current(i);
    J_n  = -Iapp./S_n;
    J_p  = Iapp./S_p;
    x_nsurf = x_navg - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
    x_psurf_set = x_pavg_set - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);

    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    eta_n    = R_g * T ./ (F * alpha_a) .* log( ( J_n + (-4*c_e*F.^2*c_nmax.^2*k_n.^2.*x_nsurf.^2 ...
       + 4*c_e*F^2*c_nmax.^2*k_n.^2.*x_nsurf+J_n.^2).^0.5 ) ./ (2*F*c_e^0.5*k_n.*(c_nmax.*x_nsurf).^0.5 .* ...
       (c_nmax-c_nmax.*x_nsurf).^0.5) ) - Iapp./S_n .* R_film;
   
    eta_p_set    = R_g * T / (F * alpha_c) .* log( ( J_p + (-4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set.^2 ...
       + 4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set+J_p.^2).^0.5 ) ./ (2*F*c_e^0.5*k_p.*(c_pmax.*x_psurf_set).^0.5 .* ...
       (c_pmax-c_pmax.*x_psurf_set).^0.5) );
   
    V_cell_set(i,:) = real(U_p_set + eta_p_set - U_n - eta_n);
    if i<(stop_cycle-29)
    weights         = 1/(0.05*sqrt(2*pi)).*exp(-(V_cell_set(i,:)-second_discharge_voltage(i)).^2/(2*0.05^2));
    weights   = weights/sum(weights) ;
    end
    x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
    x_pavg_set    = x_pavg_set - 3 * J_p / (F * R_p * c_pmax) ;
    
    if i==50
            weights_50 = weights;

    end
    if i==100
            weights_100 = weights;
    end
    if i==150
            weights_150 = weights;
    end
%     toc
end
toc
% figure
% hold on
% plot(first_discharge_time, first_discharge_voltage, 'ko');
% 
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% axis([0 4000 2.5 4])
disp('Now PF onetime calc')
k = 1;
for j = 1:50
    vec = V_cell_set(:,j);
    while vec(k,1)>2.5  && k<stop_cycle
        k=k+1;
    end
    EOD(j) = k;
    for_plot(j,1:k) = vec(1:k)';
    k=1;
end

axis([0 4300 2.5 4])

axis square
toc
[EOD ind] = sort(EOD);
weights_50s = weights_50(ind);
weights_100s = weights_100(ind);
weights_150s = weights_150(ind);

hold on
xlabel('Time(s)');
ylabel('Voltage(V)');

% plot(V_cell_set)
% weights_at_twenty = weights_at_twenty(cutoff_cycle>15); %Eliminate outliers.
% weights_at_hundred_seventy_five = weights_at_hundred_seventy_five(cutoff_cycle>15);
% weights = weights(cutoff_cycle>15);
% cutoff_cycle = cutoff_cycle(cutoff_cycle>10);
plot(second_discharge_time(EOD),weights_50s*5+2.5, 'k-', 'linewidth', 2)
plot(second_discharge_time(EOD),weights_100s*5+2.5, 'k--', 'linewidth', 2)
plot(second_discharge_time(EOD),weights_150s*5+2.5, 'k:', 'linewidth', 2)
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage, 'ko', 'linewidth', 1.5)
legend('EODV prognosis 928.6 (s)','1,853.0 (s)', '2,802.3 (s)','observations')
for j = 1:50
    plot(second_discharge_time(1:100), for_plot(j,1:100), 'color', (1-0.8*weights_50(j)/max(weights_50))*[1,1,1])
    plot(second_discharge_time(50:100), for_plot(j,50:100), 'color', (1-0.8*weights_100(j)/max(weights_100))*[1,1,1])
    plot(second_discharge_time(100:150), for_plot(j,100:150), 'color', (1-0.8*weights_150(j)/max(weights_150))*[1,1,1])
    plot(second_discharge_time(150:length(for_plot)), for_plot(j,150:end), 'color', (1-0.8*weights_150(j)/max(weights_150))*[1,1,1])
 
end
plot(second_discharge_time(EOD),weights_50s*5+2.5, 'k-', 'linewidth', 2)
plot(second_discharge_time(EOD),weights_100s*5+2.5, 'k--', 'linewidth', 2)
plot(second_discharge_time(EOD),weights_150s*5+2.5, 'k:', 'linewidth', 2)
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage, 'ko', 'linewidth', 1.5)
axis([0 4300 2.4 4])

plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
axis square
hold off



err_early = sum(weights_50.*second_discharge_time(EOD))-3653
err_late  = sum(weights_100.*second_discharge_time(EOD))-3653
err_final = sum(weights_150.*second_discharge_time(EOD))-3653
sig_early = sqrt(sum(weights_50.*(second_discharge_time(EOD) - (err_early + 3653)).^2) )
sig_late  = sqrt(sum(weights_100.*(second_discharge_time(EOD) - (err_late + 3653)).^2) )
sig_final = sqrt(sum(weights_150.*(second_discharge_time(EOD) - (err_final + 3653)).^2) )

%Now do UKF
x_navg = 0.8957971;
x_pavg = 0.5075848;
R_film = 0;
P = [0.02 0 0; 0 0.03 0; 0 0 0.01] ;
Chi = zeros(3,7);
W   = zeros(1,7);
kappa = 0.2;
S_p = 0.2536;
S_n    = 0.2447;

for i=1:100

    % temporal update requires no code.
    % measurement update
%     try
%         C = chol(P);
%     catch err
%         P = diag([0.1581 0.0001 0.01 0.033]);
%     end
    C        = chol(P);
    Chi(:,1) = [x_navg; R_film; S_p];
    Chi(:,2) = [x_navg; R_film; S_p] + sqrt(2 + kappa)*[C(1,1); 0; 0];
    Chi(:,3) = [x_navg; R_film; S_p] + sqrt(2 + kappa)*[0; C(2,2); 0];
    Chi(:,4) = [x_navg; R_film; S_p] + sqrt(2 + kappa)*[0; 0; C(3,3)];
    Chi(:,5) = [x_navg; R_film; S_p] - sqrt(2 + kappa)*[C(1,1); 0; 0];
    Chi(:,6) = [x_navg; R_film; S_p] - sqrt(2 + kappa)*[0; C(2,2); 0];
    Chi(:,7) = [x_navg; R_film; S_p] - sqrt(2 + kappa)*[0; 0; C(3,3)];   



    W(1)     = kappa / sqrt(3 + kappa);
    W(2:7)   = 1     / (2*(3+kappa));
    
    %%%%%Now the SP measurement model
    Iapp = second_discharge_current(i);
    J_n  = -Iapp/S_n;
    J_p  = Iapp./Chi(3,:);
    x_nsurf = Chi(1,:) - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
    x_psurf_set = x_pavg - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);

    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf.^0.5 - .0172./x_nsurf + ...
        .0019./x_nsurf.^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    eta_n    = R_g * T / (F * alpha_a) * log( ( J_n + (-4*c_e*F^2*c_nmax^2*k_n^2*x_nsurf.^2 ...
       + 4*c_e*F^2*c_nmax^2*k_n^2.*x_nsurf+J_n.^2).^0.5 ) ./ (2*F*c_e^0.5*k_n*(c_nmax*x_nsurf).^0.5 .* ...
       (c_nmax-c_nmax.*x_nsurf).^0.5) ) + Iapp./S_n .* Chi(2,:);
   
    eta_p_set    = R_g * T / (F * alpha_c) * log( ( J_p + (-4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set.^2 ...
       + 4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set+J_p.^2).^0.5 ) ./ (2*F*c_e^0.5*k_p.*(c_pmax.*x_psurf_set).^0.5 .* ...
       (c_pmax-c_pmax.*x_psurf_set).^0.5) );
   
    zeta = real(U_p_set + eta_p_set - U_n - eta_n);
    y_hat    = sum(W.*zeta);
    P_yy     = sum(W.*((zeta-y_hat).^2));
    y_hat_storage(i) = y_hat;
    
    P_xy     = (Chi - repmat([x_navg; R_film; S_p],1,7)) .* [W; W; W] *(zeta - y_hat)'; %P_xz is (1x1).
    K_k      = P_xy/(P_yy+0.05);
    egg = [x_navg; R_film; S_p] + K_k*(second_discharge_voltage(i) - y_hat);
    x_navg   = egg(1);
    R_film   = egg(2);
    S_p      = egg(3);
    J_p      = Iapp/S_p;
    
    err(i)   = (second_discharge_voltage(i) - y_hat);
    Q        = 0.4e-4 * ones(3,3);
    P        = P - K_k * P_xy' + Q;

    x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
    x_pavg    = x_pavg - 3 * J_p / (F * R_p * c_pmax) ;


end

x_navg
R_film
S_p




disp('Now UKF EOD calc time')
tic
zeta_storage=[];
for i=100:length(second_discharge_time)
    if x_navg <0
        x_navg = 0.01;
    end
    Iapp = second_discharge_current(i);
    J_n  = -Iapp/S_n;
    J_p  = Iapp/S_p;
    x_nsurf = x_navg - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
    x_psurf_set = x_pavg - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);

    U_n     = .7222 + .1387*x_nsurf + .029*x_nsurf^0.5 - .0172/x_nsurf + ...
        .0019/x_nsurf^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 * exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p_set     = ( -4.656 + 88.669 * x_psurf_set.^2 - 401.119 * x_psurf_set.^4 + 342.909 * ...
        x_psurf_set.^6 - 462.471 * x_psurf_set.^8 + 433.434 * x_psurf_set.^10 ) ./ ...
        ( -1 + 18.933*x_psurf_set.^2 - 79.532 * x_psurf_set.^4 + 37.311 * x_psurf_set.^6 ...
        - 73.083 * x_psurf_set.^8 + 95.96*x_psurf_set.^10 );
    eta_n    = R_g * T / (F * alpha_a) * log( ( J_n + (-4*c_e*F^2*c_nmax^2*k_n^2*x_nsurf^2 ...
       + 4*c_e*F^2*c_nmax^2*k_n^2*x_nsurf+J_n^2)^0.5 ) / (2*F*c_e^0.5*k_n*(c_nmax*x_nsurf)^0.5 * ...
       (c_nmax-c_nmax*x_nsurf)^0.5) ) + Iapp./S_n .* R_film;
   
    eta_p_set    = R_g * T / (F * alpha_c) * log( ( J_p + (-4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set.^2 ...
       + 4*c_e*F^2*c_pmax^2*k_p^2.*x_psurf_set+J_p^2).^0.5 ) ./ (2*F*c_e^0.5*k_p.*(c_pmax.*x_psurf_set).^0.5 .* ...
       (c_pmax-c_pmax.*x_psurf_set).^0.5) );
   
    zeta = real(U_p_set + eta_p_set - U_n - eta_n);
    zeta_storage(i-99) = sum(W.*zeta);
    
    x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
    x_pavg    = x_pavg - 3 * J_p / (F * R_p * c_pmax) ;
end
toc
for j = 1:length(zeta_storage)
    if zeta_storage(j) <= 2.5
        stop_here = j;
        break
    end
end
zeta_storage = zeta_storage(1:j);
figure
hold on
plot(second_discharge_time(1:100), y_hat_storage, 'k--', 'linewidth', 2)
plot(second_discharge_time(100:(100+length(zeta_storage)-1)), zeta_storage, 'k--', 'linewidth',2)
axis([0 4300 2.4 4])
xlabel('Time (s)')
ylabel('Voltage (V)')


legend('UKF prediction 1,853.0 (s)')
plot([1,4250],2.5*[1,1],'k-','linewidth',1.5) % The EODV line.
text(25,2.55,'EODV failure threshold')
plot(second_discharge_time(1:length(second_discharge_voltage)),second_discharge_voltage,'ko')
axis square
box on
error = second_discharge_time((150+length(zeta_storage)-1)) - 3653
