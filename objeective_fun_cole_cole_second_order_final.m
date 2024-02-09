function [objeective_fun_value,E_model]=objeective_fun_cole_cole_second_order_final(parameters)

%% load measured data

global freq permitivity_measured loss_factor_measured
Eo = 8.854e-12; % Permittivity of free space
%

%% Define parameter
% Einf=parameters(1);
% cond=parameters(2);

% EsCole_1=parameters(3);
% EsCole_2=parameters(4);


% tau_k_1=parameters(5);
% tau_k_2=parameters(6);


% alpha_1=parameters(7);
% alpha_2=parameters(8);

%% define model 




for f=1:1:length(freq)

tmp_1=(1j*2*pi*freq(f)*10.^parameters(5))^(1-parameters(7));
tmp_2=(1j*2*pi*freq(f)*10.^parameters(6))^(1-parameters(8));
% tmp_3=(1j*2*pi*freq(f)*10.^parameters(9))^(1-10.^parameters(13));
% tmp_4=(1j*2*pi*freq(f)*10.^parameters(10))^(1-10.^parameters(14));

E_model(f)=10.^parameters(1)+ (10.^parameters(2)/(1j*2*pi*freq(f))* Eo)+...
    +(10.^parameters(3)/(1+tmp_1))+(10.^parameters(4)/(1+tmp_2));


% loss_factor(i)=-imag(E_model)*Eo*2*pi*freq(f); % loss_factor;

end
%% for defining error and objective function 


permitivity_model=real(E_model');
loss_factor_model=imag(E_model'); % loss_factor;

permitivity_diff=sum(((permitivity_measured-permitivity_model)./median(permitivity_model)).^2);
loss_factor_diff=sum(((loss_factor_measured-loss_factor_model)./median(loss_factor_model)).^2);

 objeective_fun_value=(permitivity_diff+loss_factor_diff)/length(freq);
% objeective_fun_value=(permitivity_diff)/length(freq);
% objeective_fun_value=(loss_factor_diff)/length(freq);

% permitivity_diff=sum(abs(permitivity_measured-permitivity_model)./permitivity_measured);
% loss_factor_diff=sum(abs(loss_factor_measured-loss_factor_model)./loss_factor_measured);
% 
% objeective_fun_value=(permitivity_diff+loss_factor_diff)/length(freq);

end 