addpath('./common_function')

load network_used_in_Fig2g.mat

initial = rand(N,1)*0.8;
initial(setdiff(1:N,[1     2     3     4     6     8    10    11    12    13    14    15])) = 0;
[XX_health,X_health]=glv_Euler_time(initial,A,r,time);

X_anti = X_health;
X_anti([3     6    11    12    13]) = 0;
[XX_disease,X_disease]=glv_Euler_time(X_anti,A,r,time);

X_donor = dx_health{2}(end,:)';
X_reintro = X_disease*0.01 + X_donor*0.1;
[XX_rehealth,X_rehealth]=glv_Euler_time(X_reintro,A,r,time);

ALL = [XX_health XX_disease XX_rehealth];

MyColor = [hex2rgb('#9A6324');hex2rgb('#808000');hex2rgb('#469990')
           hex2rgb('#000075');hex2rgb('#000000');hex2rgb('#e6194B')
           hex2rgb('#f58231');hex2rgb('#ffe119');hex2rgb('#bfef45')
           hex2rgb('#3cb44b');hex2rgb('#42d4f4');hex2rgb('#4363d8')
           hex2rgb('#911eb4');hex2rgb('#f032e6');hex2rgb('#a9a9a9')];
figure;hold on;
for i = 1 : 15
    plot([time time(end)+time 2*time(end)+time],ALL(i,:),'Color',MyColor(i,:),'LineWidth',1)
end
h = plot([time time(end)+time 2*time(end)+time],ALL(4,:),'Color',[1 0 0],'LineWidth',4);

set(gca,'fontsize',14);
set(gca,'xtick',[0:30:90])
xlim([0 3*time(end)])
set(gcf,'position',[202 400 1453 322])
xlabel('time')
ylabel('Abundance')
legend(h,'C. difficile','Location','best')