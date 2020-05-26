function [A,r] = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,time,FunctionType,h1,h2,Cdiff,Cdiff_disease_abundance,Cdiff_health_abundance,select_white_black_mixed)

while 1
    [A,r] = generate_A(N,C,delta,diag,VarianceType);
    
    if strcmp(select_white_black_mixed,'mixed')
%         local_W_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
%         local_B_index = unique([Cdiff;randsample([1:1:N],randi(N))']);
        promoters  = find(A(:,Cdiff)>0);
        inhibitors = find(A(:,Cdiff)<0);
        neutrals    = find(A(:,Cdiff)==0);
        
        local_W_index = unique([Cdiff;randsample(inhibitors,randi(length(inhibitors))); randsample(promoters,randi(ceil(0.1*length(promoters)))); randsample(neutrals,randi(length(neutrals)))]);     
        local_B_index = unique([Cdiff;randsample(promoters,randi(length(promoters))); randsample(inhibitors,randi(ceil(0.1*length(inhibitors)))); randsample(neutrals,randi(length(neutrals)))]);
    else
        local_W_index = unique([Cdiff;randsample(find(A(:,Cdiff)<=0),randi(length(find(A(:,Cdiff)<=0))))]);
        local_B_index = unique([Cdiff;randsample(find(A(:,Cdiff)>=0),randi(length(find(A(:,Cdiff)>=0))))]);
    end
    initial = zeros(N,1);
    initial(local_W_index) = 0.2;
    [dx_W,dxx_W]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
    
    initial = zeros(N,1);
    initial(local_B_index) = 0.2;
    [dx_B,dxx_B]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
    
    if dxx_W(Cdiff)<Cdiff_health_abundance && dxx_B(Cdiff)>Cdiff_disease_abundance
        break;
    end
end
    
end

function [A,r] = generate_A(N,C,delta,diag,VarianceType)
r = rand(1,N);
A = zeros(N,N);
if VarianceType == 1
    zigma = 1/sqrt(N*(2+delta));
elseif VarianceType == 2
    zigma = delta;
end
for i = 1 : N
   for j = 1 : N
       if rand(1) < C
           A(i,j) =  ((zigma * randn(1)));
       end
   end
   A(i,i) = diag;
end
end