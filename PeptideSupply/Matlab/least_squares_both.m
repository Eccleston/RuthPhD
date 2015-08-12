function [err]=least_squares_both(sf, gtarget_none, gcomp_none, gtarget_ifn, gcomp_ifn, data_none, data_ifn)

%define parameters for fitting
sf1=sf(1);
sf2=sf(2);
gM1=sf(3);
gM2=sf(4);
%sf3=sf(5);
% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);

Ng1=8; Ng2=7;

gtarget_ifn=gtarget_ifn(1:7,1:7);
gcomp_ifn=gcomp_ifn(1:7,1:7);
data_ifn=data_ifn(1:7,1:7);

for i1_none=1:Ng1
    for i2_none=1:Ng1
        g1_none=gtarget_none(i1_none, i2_none);
        g2_none=gcomp_none(i1_none, i2_none);
        
        [MeP1_none(i1_none, i2_none), MeP2_none(i1_none, i2_none)]=simulateMHC(sf1*g1_none, sf2*g2_none, u1, u2, gM1);
    end
end
for i1_ifn=1:Ng2
    for i2_ifn=1:Ng2
        g1_ifn=gtarget_ifn(i1_ifn, i2_ifn);
        g2_ifn=gcomp_ifn(i1_ifn, i2_ifn);
        
        [MeP1_ifn(i1_ifn, i2_ifn), MeP2_ifn(i1_ifn, i2_ifn)]=simulateMHC_ifn(sf1*g1_ifn, sf2*g2_ifn, u1, u2, gM2);
    end
end

[p_none,e_none]=polyfit(MeP1_none, data_none,1);
err_none=e_none.normr;
%err_none = sum(((MeP1_none*sf3) - data_none).^2);
[p_ifn,e_ifn]=polyfit(MeP1_ifn, data_ifn, 1);
err_ifn=e_ifn.normr;
%err_ifn = sum(((MeP1_ifn*sf3) - data_ifn).^2);
err=err_none+err_ifn; %sum(err_none)+sum(err_ifn);

%MeP1_both=[MeP1_none(1:7,1:7), MeP1_ifn];
%data_both=[data_none(1:7,1:7), data_ifn];
%[p_both, e_both]=polyfit(MeP1_both, data_both, 1);
%err = e_both.normr;
sf
%p_none(1)=sf3;%p_both;
%p_ifn(1)=sf3;%p_both;
%p_none(2)=0;
%p_ifn(2)=0;
p_none
p_ifn
%p_both
cmap=hsv(8);
figure(8);
%subplot(1,2,1)
semilogx(gtarget_none(:,1),p_none(1)*MeP1_none(:,1)+p_none(2),gtarget_none(:,1),data_none(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_none(:,2),p_none(1)*MeP1_none(:,2)+p_none(2),gtarget_none(:,2),data_none(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,3),p_none(1)*MeP1_none(:,3)+p_none(2),gtarget_none(:,3),data_none(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,4),p_none(1)*MeP1_none(:,4)+p_none(2),gtarget_none(:,4),data_none(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,5),p_none(1)*MeP1_none(:,5)+p_none(2),gtarget_none(:,5),data_none(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,6),p_none(1)*MeP1_none(:,6)+p_none(2),gtarget_none(:,6),data_none(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,7),p_none(1)*MeP1_none(:,7)+p_none(2),gtarget_none(:,7),data_none(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none(:,8),p_none(1)*MeP1_none(:,8)+p_none(2),gtarget_none(:,8),data_none(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
%legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
figure(9);
%subplot(1,2,2)
semilogx(gtarget_ifn(:,1),p_ifn(1)*MeP1_ifn(:,1)+p_ifn(2),gtarget_ifn(:,1),data_ifn(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_ifn(:,2),p_ifn(1)*MeP1_ifn(:,2)+p_ifn(2),gtarget_ifn(:,2),data_ifn(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn(:,3),p_ifn(1)*MeP1_ifn(:,3)+p_ifn(2),gtarget_ifn(:,3),data_ifn(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn(:,4),p_ifn(1)*MeP1_ifn(:,4)+p_ifn(2),gtarget_ifn(:,4),data_ifn(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn(:,5),p_ifn(1)*MeP1_ifn(:,5)+p_ifn(2),gtarget_ifn(:,5),data_ifn(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn(:,6),p_ifn(1)*MeP1_ifn(:,6)+p_ifn(2),gtarget_ifn(:,6),data_ifn(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn(:,7),p_ifn(1)*MeP1_ifn(:,7)+p_ifn(2),gtarget_ifn(:,7),data_ifn(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
%semilogx(gtarget(:,8),p(1)*MeP1(:,8)+p(2),gtarget(:,8),data(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' '})
set(gca, 'FontSize',16);
hold off
return

