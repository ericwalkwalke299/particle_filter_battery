function [states,n_iter] = my_lsqnonlin_He(obs)
% This function performs non-linear least squares with 
% He, et al's model.
IG  = [-9.86e-7;5.752e-2;8.983e-1;-8.34e-4];   
obs = obs';
u   = (1:length(obs))';
n_iter=0;
delta_IG=IG;
while n_iter < 200 && sum(abs(delta_IG./IG) > 1e-4)/4
    jac = [exp(IG(2)*u), IG(1)*IG(2)*exp(IG(2)*u),...
    exp(IG(4)*u), IG(3) * IG(4) * exp(IG(4)*u)];
    obs_model = IG(1)*exp(IG(2)*u) + IG(3) * exp(IG(4)*u);
    delta_IG = (jac'*jac)\jac'*(obs-obs_model);
    IG = IG + delta_IG; n_iter=n_iter+1;
end

states = IG;