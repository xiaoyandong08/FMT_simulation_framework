function FMT_simulation_Framework
addpath('./common_function')

N = 100; % community size, # of species
C = 0.4; % probability of effect from species-i to species-j
delta = 0.2; % variance of a_ij 
diag = -1.0; % a_ii
VarianceType = 2; % VarianceType = 1;a_ij \sim \mathcal N(0,1/sqrt(N*(2+delta))^); VarianceType = 2; a_ij \sim \mathcal N(0,delta^2);
time = [0:0.1:30];
FunctionType = 1; % FunctionType = 1, GLV, =2, Holling Type II; = 3, DeAngelis-Beddington; =4, Crowley-Martin
h1 = 0.1; % parameter in Holling Type II, DeAngelis-Beddington, and Crowley-Martin
h2 = 0.1; % parameter in Holling Type II, DeAngelis-Beddington, and Crowley-Martin

Cdiff_disease_abundance = 0.5;
Cdiff_health_abundance = 1e-4;
select_white_black_mixed = 'mixed';
Cdiff = 1;

Disease_threshold = Cdiff_disease_abundance;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate A matrix and growth rate for a microbial community
[A,r] = Generate_Network_A_of_BandW(N,C,delta,diag,VarianceType,time,FunctionType,h1,h2,Cdiff,Cdiff_disease_abundance,Cdiff_health_abundance,select_white_black_mixed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate X_disease sample

min_threshold_rCDI = 10;
max_threshold_rCDI = 15;
[XX_disease,X_disease,XX_health,X_health] = Generate_disease_sample(A,r,time,FunctionType,h1,h2,min_threshold_rCDI,max_threshold_rCDI,Disease_threshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate X_donor sample
min_donor_species = 60;
max_donor_species = 100;

[XX_donor,X_donor] = Generate_donor_samples(N,A,r,Cdiff,Cdiff_health_abundance,time,FunctionType,h1,h2,min_donor_species,max_donor_species);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transplanting process

X_reintro_donor = X_disease + X_donor;
[XX_FMT,X_FMT]=glv_Euler_type(X_reintro_donor,A,r,time,FunctionType,h1,h2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL = [XX_health' XX_disease' XX_FMT'];

MyColor = lines(N);
figure;hold on;
for i = 1 : N
    plot([time time(end)+time 2*time(end)+time],ALL(i,:),'Color',MyColor(i,:),'LineWidth',1)
end
h = plot([time time(end)+time 2*time(end)+time],ALL(Cdiff,:),'Color',[1 0 0],'LineWidth',4);

set(gca,'fontsize',14);
set(gca,'xtick',[0:30:90])
xlim([0 3*time(end)])
set(gcf,'position',[202 400 1453 322])
xlabel('time')
ylabel('Abundance')
legend(h,'C. difficile','Location','best')
end

