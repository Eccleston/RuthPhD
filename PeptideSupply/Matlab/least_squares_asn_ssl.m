function [err]=least_squares_asn_ssl(sf, gtarget, gcomp,data, gtarget_2, gcomp_2, data_2)

%define parameters for fitting
sf1=sf(1);
sf2=sf(2);
%gM1=sf(3);
%upreg1=sf(4);
%upreg2=sf(5);
%gM2=sf(4);
%sf3=sf(5);
% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);

Ng1=8; Ng2=8;

%gtarget_ifn=gtarget_ifn(1:7,1:7);
%gcomp_ifn=gcomp_ifn(1:7,1:7);
%data_ifn=data_ifn(1:7,1:7);

for i1_none=1:Ng1
    for i2_none=1:Ng1
        g1_none=gtarget(i1_none, i2_none);
        g2_none=gcomp(i1_none, i2_none);
        
        [MeP1_none(i1_none, i2_none), MeP2_none(i1_none, i2_none)]=simulateMHC_asn_ssl(sf1*g1_none, sf2*g2_none, u1, u2, 105.5);
    end
end
for i1_none_2=1:Ng2
    for i2_none_2=1:Ng2
        g1_none_2=gtarget_2(i1_none_2, i2_none_2);
        g2_none_2=gcomp_2(i1_none_2, i2_none_2);
        
        [MeP1_none_2(i1_none_2, i2_none_2), MeP2_none_2(i1_none_2, i2_none_2)]=simulateMHC_asn_ssl(sf1*g1_none_2, sf2*g2_none_2, u1, u2, 105.5);
    end
end

validMeP1_none = ~isnan(MeP1_none);
validdata = ~isnan(data);
validdataBoth = validMeP1_none & validdata;
MeP1_none_valid = MeP1_none(validdataBoth) ;
data_valid = data(validdataBoth);

validMeP2_none_2 = ~isnan(MeP2_none_2);
validdata_2 = ~isnan(data_2);
validdataBoth_2 = validMeP2_none_2 & validdata_2;
MeP2_none_2_valid = MeP2_none_2(validdataBoth_2) 
data_valid_2 = data_2(validdataBoth_2);


[p_MeP1,e_MeP1]=polyfit(MeP1_none_valid, data_valid,1);
err_MeP1=e_MeP1.normr;

%err_none = sum(((MeP1_none*sf3) - data).^2);
[p_MeP2,e_MeP2]=polyfit(MeP2_none_2_valid, data_valid_2, 1);
err_MeP2=e_MeP2.normr;
%err_ifn = sum(((MeP1_ifn*sf3) - data_ifn).^2);
err=err_MeP1+err_MeP2; %sum(err_none)+sum(err_ifn);

%MeP1_both=[MeP1_none(1:7,1:7), MeP1_none_2];
%data_both=[data(1:7,1:7), data_ifn];
%[p_both, e_both]=polyfit(MeP1_both, data_both, 1);
%err = e_both.normr;
sf
%p_MeP1=p_both;
%p_ifn=p_both;
%p_MeP1(2)=0;
%p_ifn(2)=0;
%p_MeP1
%p_ifn
p_MeP1
p_MeP2
cmap=hsv(8);
figure(8);
%subplot(1,2,1)
semilogx(gtarget(:,1),p_MeP1(1)*MeP1_none(:,1)+p_MeP1(2),gtarget(:,1),data(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget(:,2),p_MeP1(1)*MeP1_none(:,2)+p_MeP1(2),gtarget(:,2),data(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,3),p_MeP1(1)*MeP1_none(:,3)+p_MeP1(2),gtarget(:,3),data(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,4),p_MeP1(1)*MeP1_none(:,4)+p_MeP1(2),gtarget(:,4),data(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,5),p_MeP1(1)*MeP1_none(:,5)+p_MeP1(2),gtarget(:,5),data(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,6),p_MeP1(1)*MeP1_none(:,6)+p_MeP1(2),gtarget(:,6),data(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,7),p_MeP1(1)*MeP1_none(:,7)+p_MeP1(2),gtarget(:,7),data(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget(:,8),p_MeP1(1)*MeP1_none(:,8)+p_MeP1(2),gtarget(:,8),data(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
figure(9);
%subplot(1,2,2)
semilogx(gtarget_2(:,1),p_MeP2(1)*MeP2_none_2(:,1)+p_MeP2(2),gtarget_2(:,1),data_2(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_2(:,2),p_MeP2(1)*MeP2_none_2(:,2)+p_MeP2(2),gtarget_2(:,2),data_2(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,3),p_MeP2(1)*MeP2_none_2(:,3)+p_MeP2(2),gtarget_2(:,3),data_2(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,4),p_MeP2(1)*MeP2_none_2(:,4)+p_MeP2(2),gtarget_2(:,4),data_2(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,5),p_MeP2(1)*MeP2_none_2(:,5)+p_MeP2(2),gtarget_2(:,5),data_2(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,6),p_MeP2(1)*MeP2_none_2(:,6)+p_MeP2(2),gtarget_2(:,6),data_2(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,7),p_MeP2(1)*MeP2_none_2(:,7)+p_MeP2(2),gtarget_2(:,7),data_2(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_2(:,8),p_MeP2(1)*MeP2_none_2(:,8)+p_MeP2(2),gtarget_2(:,8),data_2(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
return

