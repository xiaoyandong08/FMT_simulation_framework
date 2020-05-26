function [XX_donor,X_donor] = Generate_donor_samples(N,A,r,disordered_species,X_health_target,time,FunctionType,h1,h2,min_donor_species,max_donor_species)
Cdiff = 1;
promoters  = find(A(:,Cdiff)>0);
inhibitors = find(A(:,Cdiff)<0);
neutrals    = find(A(:,Cdiff)==0);

flag = 1;
while flag == 1
    index = unique([Cdiff randsample(inhibitors,randi(length(inhibitors)))' randsample(neutrals,randi(length(neutrals)))' randsample(promoters,randi(length(promoters)))']);

    donor_species_index = zeros(1,N);
    donor_species_index(index) = index;
    donor_species_index(disordered_species) = disordered_species;
    
    initial = zeros(N,1);
    initial(donor_species_index~=0) = 0.2;
    [XX_donor,X_donor]=glv_Euler_type(initial,A,r,time,FunctionType,h1,h2);
    donor_species_richness = sum(X_donor>0);

    if X_donor(disordered_species)<X_health_target && donor_species_richness>min_donor_species && donor_species_richness<max_donor_species
        flag = 0;
    end
end

end

