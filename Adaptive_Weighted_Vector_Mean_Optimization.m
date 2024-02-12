%% Adaptive Weighted Vector Mean Optimization 

% This code is inspired by the following work and needes to be cited their works. 
% Iman Ahmadianfar, Ali asghar Heidari, Saeed Noushadian, Huiling Chen, Amir H. Gandomi  
% INFO: An Efficient Optimization Algorithm based on Weighted Mean of Vectors
% Expert Systems With Applications, 116516, 2022, doi: https://doi.org/10.1016/j.eswa.2022.116516
% http://www.aliasgharheidari.com/INFO.html


%
function [Best_Cost,Best_X,Convergence_curve]=Adaptive_Weighted_Vector_Mean_Optimization(nP,MaxIt,lb,ub,dim,fobj)
%% Adaptive_initilization

        Cost=zeros(nP,1);
        M=zeros(nP,1);
        
        X= adaptive_initilization(nP,dim,lb,ub,fobj);  %% using Adaptive initialization 
        
        for i=1:nP
           Cost(i) = fobj(X(i,:)); 
           M(i)=Cost(i);
        end
        
        [~, ind]=sort(Cost);
        Best_X = X(ind(1),:);
        Best_Cost = Cost(ind(1));
        
        Worst_Cost = Cost(ind(end));
        Worst_X = X(ind(end),:);
        
        I=randi([2 5]);
        Better_X=X(ind(I),:);
        Better_Cost=Cost(ind(I));
        Max_iter_diversity_factor=floor(MaxIt*0.2);
        
%% Main Loop of Adaptive_Weighted_Vector_Mean_Optimization
        for it=1:MaxIt
            alpha=2*exp(-4*(it/MaxIt));                                                          
            
            M_Best=Best_Cost;
            M_Better=Better_Cost;
            M_Worst=Worst_Cost;
            
            for i=1:nP
                
               % Updating rule stage
                del=2*rand*alpha-alpha;                                                         
                sigm=2*rand*alpha-alpha;                                                                                                
                                                                 
                % Select three random solution
                A1=randperm(nP);
                A1(A1==i)=[];
                a=A1(1);b=A1(2);c=A1(3);
                
                e=1e-25;
                epsi=e*rand;
                
                omg = max([M(a) M(b) M(c)]);
                MM = [(M(a)-M(b)) (M(a)-M(c)) (M(b)-M(c))];
                
                W(1) = cos(MM(1)+pi)*exp(-abs(MM(1)/omg));                                          
                W(2) = cos(MM(2)+pi)*exp(-abs(MM(2)/omg));                                         
                W(3)= cos(MM(3)+pi)*exp(-abs(MM(3)/omg));                                            
                Wt = sum(W);
                
                WM1 = del.*(W(1).*(X(a,:)-X(b,:))+W(2).*(X(a,:)-X(c,:))+ ...                      
                    W(3).*(X(b,:)-X(c,:)))/(Wt+1)+epsi;
                
                omg = max([M_Best M_Better M_Worst]);
                MM = [(M_Best-M_Better) (M_Best-M_Better) (M_Better-M_Worst)];
                
                W(1) = cos(MM(1)+pi)*exp(-abs(MM(1)/omg));                                        
                W(2) = cos(MM(2)+pi)*exp(-abs(MM(2)/omg));                                           
                W(3) = cos(MM(3)+pi)*exp(-abs(MM(3)/omg));                                           
                Wt = sum(W);
                
                WM2 = del.*(W(1).*(Best_X-Better_X)+W(2).*(Best_X-Worst_X)+ ...                
                    W(3).*(Better_X-Worst_X))/(Wt+1)+epsi;
                
                % Determine MeanRule 
                r = unifrnd(0.1,0.5);
                MeanRule = r.*WM1+(1-r).*WM2;                                                    
                
                if rand<0.5
                    z1 = X(i,:)+sigm.*(rand.*MeanRule)+randn.*(Best_X-X(a,:))/(M_Best-M(a)+1);
                    z2 = Best_X+sigm.*(rand.*MeanRule)+randn.*(X(a,:)-X(b,:))/(M(a)-M(b)+1);
                else                                                                              
                    z1 = X(a,:)+sigm.*(rand.*MeanRule)+randn.*(X(b,:)-X(c,:))/(M(b)-M(c)+1);
                    z2 = Better_X+sigm.*(rand.*MeanRule)+randn.*(X(a,:)-X(b,:))/(M(a)-M(b)+1);
                end
                
               % Vector combining stage
                u=zeros(1,dim);
                for j=1:dim
                    mu = 0.05*randn;
                    if rand <0.5 
                        if rand<0.5
                            u(j) = z1(j) + mu*abs(z1(j)-z2(j));                                  
                        else
                            u(j) = z2(j) + mu*abs(z1(j)-z2(j));                                  
                        end
                    else
                        u(j) = X(i,j);                                                            
                    end
                end
                
                % Local search stage
                if rand<0.5
                    L=rand<0.5;v1=(1-L)*2*(rand)+L;v2=rand.*L+(1-L);                             
                    Xavg=(X(a,:)+X(b,:)+X(c,:))/3;                                               
                    phi=rand;
                    Xrnd = phi.*(Xavg)+(1-phi)*(phi.*Better_X+(1-phi).*Best_X);                  
                    Randn = L.*randn(1,dim)+(1-L).*randn;
                    if rand<0.5
                        u = Best_X + Randn.*(MeanRule+randn.*(Best_X-X(a,:)));                    
                    else
                        u = Xrnd + Randn.*(MeanRule+randn.*(v1*Best_X-v2*Xrnd));                 
                    end
                    
                end
                
                % Check if new solution go outside the search space and bring them back
                New_X= BC(u,lb,ub);
                New_Cost = fobj(New_X);
                
                if New_Cost<Cost(i)
                    X(i,:)=New_X;
                    Cost(i)=New_Cost;
                    M(i)=Cost(i);
                    if Cost(i)<Best_Cost
                        Best_X=X(i,:);
                        Best_Cost = Cost(i);
                    end
                end
            end


             %% Reintroduce diversity by reinitializing a portion of the population using Adaptive initialization 
       
        if mod(it,Max_iter_diversity_factor) == 0

            portion_to_reinitialize = floor(0.2 * nP);  %# reinitialize 20% of the population
            tmp=adaptive_initilization(portion_to_reinitialize,dim,lb,ub,fobj);  %% using Adaptive initialization 
            % X[(nP-portion_to_reinitialize):nP, :] = initialize_population_MH(portion_to_reinitialize, dim, ub, lb, fobj);
            X(nP-portion_to_reinitialize+1:nP,:)=tmp;

        end
          

        
            % Determine the worst solution
            [~, ind]=sort(Cost);
            Worst_X=X(ind(end),:);
            Worst_Cost=Cost(ind(end));
            % Determine the better solution
            I=randi([2 5]);
            Better_X=X(ind(I),:);
            Better_Cost=Cost(ind(I));

            % Update Convergence_curve
            Convergence_curve(it)=Best_Cost;

            % Show Iteration Information
%              disp(['Iteration ' num2str(it) ',: Best Cost = ' num2str(Best_Cost)]);
        end

end


function X = BC(X,lb,ub) 
Flag4ub=X>ub;
Flag4lb=X<lb;
X=(X.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end    
    
    
    
    
    
    




