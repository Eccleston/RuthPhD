function [val]=least_squares_fun(sf,gtarget,gcomp,data)

sf1=sf(1);
sf2=sf(2);
sf3=sf(3);
tf=sf(4);
% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);

Ng=8;


for i1 = 1:Ng
    for i2 = 1:Ng
  g1 = gtarget(i1,i2); % 
   g2 = gcomp(i1,i2);
    [MeP1(i1,i2),MeP2(i1,i2)] = simulateMHC(sf1*g1,sf2*g2,u1,u2,tf);
    %resid=sum((sf2*MeP1(:,i1)-data(:,i1)).^2);
  end
end
%for i=1:8
resid =sum((sf3*MeP1(:)-data(:)).^2);
%end
val= nansum(resid);
%sf
cmap=hsv(8);
figure(5);
semilogx(gtarget(:,1),sf3*MeP1(:,1),gtarget(:,1),data(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget(:,2),sf3*MeP1(:,2),gtarget(:,2),data(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,3),sf3*MeP1(:,3),gtarget(:,3),data(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,4),sf3*MeP1(:,4),gtarget(:,4),data(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,5),sf3*MeP1(:,5),gtarget(:,5),data(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,6),sf3*MeP1(:,6),gtarget(:,6),data(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,7),sf3*MeP1(:,7),gtarget(:,7),data(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,8),sf3*MeP1(:,8),gtarget(:,8),data(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
set(gca, 'FontSize',16);
hold off
end