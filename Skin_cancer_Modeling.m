
%% paper Ref: Awal et al, Adaptive Weighted Vector Means Optimization for Healthy and Malignant Skin  Modeling at Microwave Frequencies Using Clinical Data
%IEEE JOURNAL OF ELECTROMAGNETICS, RF, AND MICROWAVES IN MEDICINE AND
%BIOLOGY, 2024

clear all
 close all
% clc
rng ("default")
global freq permitivity_measured loss_factor_measured  
Eo = 8.854e-12; % Permittivity of free space


%% load measured data
  load('Skin_dielectric_data.mat')


% take 0.3ghz to 14 ghz 




nP=500;          % Number of Population


MaxIt=100;      % Maximum number of iterations; in paper it is 300. 

%
%% Define parameter

Einf=[0.001; 5] ;
cond=[  -3 ;  0.1 ] ; 

EsCole_1=[ 0.5; 3] ; 
EsCole_2=[ 0.5 ; 3] ;


 tau_k_1=[  -14;  -8] ; 
 tau_k_2=[ -14; -8] ; 


alpha_1=[0.0; 1] ; 

 alpha_2= [ 0.0; 1]; 

LBUB=[Einf  cond   EsCole_1 EsCole_2 tau_k_1 tau_k_2 alpha_1 alpha_2];





lb=LBUB(1,:);

ub=LBUB(2,:);
dim=length(lb);


tStart = cputime;
[Best_fitness,BestPositions,Convergence_curve] = Adaptive_Weighted_Vector_Mean_Optimization(nP,MaxIt,lb,ub,dim,@objeective_fun_cole_cole_second_order_final);
tEnd = cputime - tStart

BestPositions
Best_fitness


Best_param(1:6)=10.^(BestPositions(1:6));
Best_param(7:8)=BestPositions(7:8);
Best_param
%% Draw objective space



figure,
hold on
semilogy(Convergence_curve,'Color','r','LineWidth',4);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('INFO')



[objeective_fun_value,E_model]=objeective_fun_cole_cole_second_order_final(BestPositions);

permitivity_model=real(E_model);
loss_factor_model=-imag(E_model);

figure


plot(freq./1e9,permitivity_measured,'ob')
hold on
plot(freq./1e9,permitivity_model,'-r','LineWidth',2)

legend('Measured data','Model')
xlabel('Frequency (GHz)')
ylabel('Dielectric Constant')


figure

plot(freq./1e9,loss_factor_measured,'ob')
hold on
plot(freq./1e9,loss_factor_model,'-r','LineWidth',2)
legend('Measured data','Model')
xlabel('Frequency (GHz)')
ylabel('Loss factor')
