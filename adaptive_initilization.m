
% This function initialize the first population of search agents
% function X=initialization(nP,dim,ub,lb)
% 
% Boundary_no= size(ub,2); % numnber of boundaries
% 
% % If the boundaries of all variables are equal and user enter a signle
% % number for both ub and lb
% if Boundary_no==1
%     X=rand(nP,dim).*(ub-lb)+lb;
% end
% 
% % If each variable has a different lb and ub
% if Boundary_no>1
%     for i=1:dim
%         ub_i=ub(i);
%         lb_i=lb(i);
%         X(:,i)=rand(nP,1).*(ub_i-lb_i)+lb_i;
%     end
% end
function Xint_Adaptive= adaptive_initilization(N,D,lb,ub,fobj)

global freq permitivity_measured loss_factor_measured

   X=initialization(2*N,D,ub,lb);
    Xint=zeros(N,D);
    Xint_Adaptive=zeros(N,D);
    [rw,cl]=size(Xint);
     if length(lb)  ==  1  
         lb  =  repmat( lb ,  [rw,cl] ) ; 
         ub  =  repmat( ub ,  [rw,cl]) ; 
     end


        
    for j=1:cl
    
        temp=linspace(lb(j),ub(j),N);
        Xint(:,j)=temp(:);
        
    end

    population = initialize_population_DRAM(N,D, lb, ub,fobj);
    Xint_final=[X;Xint;population];

    [r,c]=size(Xint_final);
        
    for i=1:r
        
        Cost(i) = fobj(Xint_final(i,:)); 
    end
          
           
        [~, ind]=sort(Cost);

         Xint_Adaptive = Xint_final(ind(1:N),:);

end


function X=initialization(nP,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(nP,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(nP,1).*(ub_i-lb_i)+lb_i;
    end
end

end






% 
% end 

% x=linspace(-1,1,10);