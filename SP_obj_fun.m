function obj = SP_obj_fun(pars,time,current,voltage)

S_n    = 0.2607;
S_p    = 0.2571;
k_n    = 37.4312e-12;
k_p    = 17.4733e-12;
R_n    = 2e-6;
R_p    = 2e-6;
D_n    = 29.0798e-15;
D_p    = 27.9034e-15;
c_nmax = 30074.5;
c_pmax = 51563.5;
c_e    = 1000;
x_navg = 0.9388;
x_pavg = 0.5171;
T      = 298.15;
R_g    = 8.3143;
F      = 96487;
alpha_a  = 0.5;
alpha_c  = 0.5;

if length(pars) > 1
    S_n = pars(1);
    S_p = pars(2);
end
if length(pars) > 2
    x_navg = pars(3);
    D_pavg = pars(4);
end
if length(pars) > 4
    k_n = pars(5);
    k_p = pars(6);
end
if length(pars) > 6
    x_navg = pars(7);
    x_pavg = pars(8);
end
if length(pars) > 8
    c_nmax = pars(9);
    c_pmax = pars(10);
end
if length(pars) > 10
    alpha_a = 0.5;
    alpha_c = 0.5;
end
if length(pars) > 12
    R_n = pars(13);
    R_p = pars(14);
end
if length(pars) > 14
    c_e = pars(15);
    T   = pars(16);
end
V_cell = [];  % The model returns voltage, which is diplayed
% in plots.  The vector is initialized before assigning entries.  
% Likewise, initialize the vectors for states of charge.
x_navg_vec = [];
x_pavg_vec = [];

% The model is put into motion.  
for i = 1:length(time)  % The model is going to calculate, 
    % one time point at an iteration, forward to the end.
    % Calculate the voltage at the current time point, first.  
    % Then, thedynamic model will reach ahead to prepare the 
    % changing states of charge for the next loop, the next 
    % time point. 
    
    Iapp    = current(i); % Assign the current at the present
    % time point, so it's less bulky in the equations.
    J_n     = -Iapp / S_n;
    J_p     = Iapp / S_p;
    x_nsurf = x_navg - ( J_n * R_n ) / ( 5 * F * D_n * c_nmax);
    x_psurf = x_pavg - ( J_p * R_p ) / ( 5 * F * D_p * c_pmax);
    
    % Now, we have enough for the open circuit potentials.
    U_n     = .8214 + .1387*x_nsurf + .029*x_nsurf^0.5 - .0172/...
        x_nsurf + ...
        .0019/x_nsurf^1.5 + .2808 * exp(0.9-15*x_nsurf) -.7984 *...
        exp (...
        0.4465*x_nsurf - 0.4108);
    
    U_p     = ( -4.8811 + 88.669 * x_psurf^2 - 401.119 * x_psurf^4 +...
        342.909 * ...
        x_psurf^6 - 462.471 * x_psurf^8 + 433.434 * x_psurf^10 ) / ...
        ( -1 + 18.933*x_psurf^2 - 79.532 * x_psurf^4 + 37.311 *...
        x_psurf^6 ...
        - 73.083 * x_psurf^8 + 95.96*x_psurf^10 );
    
   
   eta_n    = R_g * T / (F * alpha_a) * log( ( J_n + (-4*c_e*F^2*...
       c_nmax^2*k_n^2*x_nsurf^2 ...
       + 4*c_e*F^2*c_nmax^2*k_n^2*x_nsurf+J_n^2)^0.5 ) / (2*F*...
       c_e^0.5*k_n*(c_nmax*x_nsurf)^0.5 * ...
       (c_nmax-c_nmax*x_nsurf)^0.5) );
   
   eta_p    = R_g * T / (F * alpha_c) * log( ( J_p + (-4*c_e*...
       F^2*c_pmax^2*k_p^2*x_psurf^2 ...
       + 4*c_e*F^2*c_pmax^2*k_p^2*x_psurf+J_p^2)^0.5 ) / (2*...
       F*c_e^0.5*k_p*(c_pmax*x_psurf)^0.5 * ...
       (c_pmax-c_pmax*x_psurf)^0.5) );
   
   % Now, the model returns its voltage.
   V_cell(i) = U_p + eta_p - U_n - eta_n;
   
   % Prepare the state of charge for the next iteration,
   %  based upon the present current and the time step to come.
   
   if i<length(time)  % The conditional statement is 
       % necessary because at the very end, 'i + 1' is out of bounds
       % of the data vector.
       t_step    = time(i+1) - time(i);
   end % t_step will be left as the last time step, when the end 
   % of the data vector has passed.
   
   % Before changing the SOC, save the current point, for plotting.
   x_navg_vec(i) = x_navg;
   x_pavg_vec(i) = x_pavg;
   
   x_navg    = x_navg - 3 * J_n / (F * R_n * c_nmax) ;
   x_pavg    = x_pavg - 3 * J_p / (F * R_p * c_pmax) ;
   
end

obj = voltage - V_cell; % Change to (1:50) for 928.6s estimate.