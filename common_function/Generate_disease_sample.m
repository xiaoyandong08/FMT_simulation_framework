function [XX_disease,X_disease,XX_health,X_health] = Generate_disease_sample(A,r,time,FunctionType,h1,h2,min_threshold,max_threshold,Disease_threshold)
Cdiff = 1;
N = size(A,1);

promoters  = find(A(:,Cdiff)>0);
inhibitors = find(A(:,Cdiff)<0);
neutrals    = find(A(:,Cdiff)==0);
disordered_species = Cdiff;

flag = 1;
while flag == 1
    count = 0;
    [XX_health,X_health] = Generate_donor_samples(N,A,r,disordered_species,1e-4,time,FunctionType,h1,h2,50,100);
%     while 1
%         index = unique([Cdiff randsample(inhibitors,randi(length(inhibitors)))' randsample(neutrals,randi(length(neutrals)))' randsample(promoters,randi(length(promoters)))']);
%         patient_species_index = zeros(1,N);
%         patient_species_index(index) = index;
% 
%         initial = zeros(N,1);
%         initial(patient_species_index~=0) = 0.2;
%         [XX_health,X_health]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
% 
%         if X_health(disordered_species)<1e-4
%             break;
%         end
%     end
    while 1
        
        Delete =setdiff([randsample(inhibitors,randi(length(inhibitors)))' randsample(neutrals,randi(length(neutrals)))'],Cdiff);        
        X_anti = X_health;
        X_anti(Delete) = 0;
        [XX_disease,X_disease]=glv_Euler_type(X_anti,A,r,time,FunctionType,h1,h2);
        Patient_species_richness = sum(X_disease~=0);
        count = count + 1;
        if X_disease(disordered_species) - X_health(disordered_species) > Disease_threshold && Patient_species_richness>min_threshold && Patient_species_richness<max_threshold
            flag = 0;
            break;
        end
        if count>200
            break;
        end
    end
end
end
