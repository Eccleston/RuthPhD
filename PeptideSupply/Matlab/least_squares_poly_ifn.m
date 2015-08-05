function [err]=least_squares_poly_ifn(sf,gtarget,gcomp,data)

%sf1=sf(1);
%sf2=sf(2);

upreg=sf(1);

% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);
gtarget=gtarget(1:7,1:7);
gcomp=gcomp(1:7,1:7);
data=data(1:7,1:7);
Ng=7;
%Ng2=8;
for i1 = 1:Ng
    for i2 = 1:Ng
  g1 = gtarget(i1,i2); % 
   g2 = gcomp(i1,i2);
    % [MeP1(i1,i2),MeP2(i1,i2)] = simulateMHC_ifn(sf1*g1,sf2*g2,u1,u2,upreg);
     [MeP1(i1,i2),MeP2(i1,i2)] = simulateMHC_ifn(g1,g2,u1,u2,upreg);
    %resid=sum((sf2*MeP1(:,i1)-data(:,i1)).^2);
  end
end

%for i=1:8
%[p(i,:),e(i)]=polyfit(MeP1(:,i),data(:,i),1);
%resid(i) =nansum(((p(i,1)*MeP1(:,i)+p(i,2))-data(:,i)).^2);
%MeTarget(:,i)= p(i,1).*MeP1(:,i)+p(i,2);
%end
%MeP1=MeP1*upreg;
%validdata1 = ~isnan(MeP1);
%validdata2 = ~isnan(data)
%validdataBoth = validdata1 & validdata2
%MeP1valid = MeP1(validdataBoth) 
%datavalid = data(validdataBoth) 

[p,e] = polyfit(MeP1,data,1);
err = e.normr;
%err=nansum(err);
%sf3=p(1);
%MeP1=MeP1*upreg;
%p=[0.0011 361.8626];%[0.0048, 338.7170];%[0.0015 343.9527];
%p=p*upreg;
%resid =nansum(((p(1)*MeP1(:)+p(2))-data(:)).^2);% + nansum(((p(1)*MeP1(1:6,2)+p(2))-data(1:6,2)).^2) +nansum(((p(1)*MeP1(:,3:8)+p(2))-data(:,3:8)).^2);%+sum((p(1)*MeP1(:,3:8)+p(2))-data(:,3:8)).^2
%end
%err= sum(resid);
sf
p1=p(1)
p2=p(2)
cmap=hsv(8);
figure(9);
semilogx(gtarget(:,1),p(1)*MeP1(:,1)+p(2),gtarget(:,1),data(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget(:,2),p(1)*MeP1(:,2)+p(2),gtarget(:,2),data(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,3),p(1)*MeP1(:,3)+p(2),gtarget(:,3),data(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,4),p(1)*MeP1(:,4)+p(2),gtarget(:,4),data(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,5),p(1)*MeP1(:,5)+p(2),gtarget(:,5),data(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,6),p(1)*MeP1(:,6)+p(2),gtarget(:,6),data(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,7),p(1)*MeP1(:,7)+p(2),gtarget(:,7),data(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
%semilogx(gtarget(:,8),p(1)*MeP1(:,8)+p(2),gtarget(:,8),data(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
%legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' '})

set(gca, 'FontSize',16);
hold off
%plot(gtarget,p(1)*MeP1+p(2),'r--')
%%hold on
%plot(gtarget,data,'kx')
%hold off
return