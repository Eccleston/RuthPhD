function [err]=least_squares_asn_ssl_ifn_check(sf,gtarget_none_1,gcomp_none_1,data_none_1,gtarget_none_2,gcomp_none_2, data_none_2,gtarget_ifn1_1,gcomp_ifn1_1,data_ifn1_1,gtarget_ifn1_2,gcomp_ifn1_2, data_ifn1_2)

%define parameters for fitting
sf1=sf(1);
sf2=sf(2);
sf3=sf(3);
upfactor=sf(4);
gM=150.5;
% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);

Ng1=8; Ng2=8;


for i1=1:Ng1
    for i2=1:Ng1
        %g1_none=gtarget(i1, i2);
        %g2_none=gcomp(i1, i2);
        
        [MeP1_none_1(i1, i2), MeP2_none_1(i1, i2)]=simulateMHC_asn_ssl_check(sf1*gtarget_none_1(i1, i2), sf2*gcomp_none_1(i1, i2), u1, u2, gM,0,1);
        [MeP1_none_2(i1, i2), MeP2_none_2(i1, i2)]=simulateMHC_asn_ssl_check(sf1*gtarget_none_2(i1, i2), sf2*gcomp_none_2(i1, i2), u1, u2, gM,0,1);
        [MeP1_ifn1_1(i1, i2), MeP2_ifn1_1(i1, i2)]=simulateMHC_asn_ssl_check(sf1*gtarget_ifn1_1(i1, i2), sf2*gcomp_ifn1_1(i1, i2), u1, u2, gM, 1,upfactor);
        [MeP1_ifn1_2(i1, i2), MeP2_ifn1_2(i1, i2)]=simulateMHC_asn_ssl_check(sf1*gtarget_ifn1_2(i1, i2), sf2*gcomp_ifn1_2(i1, i2), u1, u2, gM, 1,upfactor);
    end
end


%validMeP1_none_1 = ~isnan(MeP1_none_1);
%validdata_none_1 = ~isnan(data_none_1);
%validdataBoth_none_1 = validMeP1_none_1 & validdata_none_1;
%MeP1_none_valid_1 = MeP1_none_1(validdataBoth_none_1) ;
%data_valid_none_1 = data_none_1(validdataBoth_none_1);

%[p_MeP1_none_1,e_MeP1_none_1]=polyfit(MeP1_none_valid_1, data_valid_none_1,1);
%err_MeP1_none_1=e_MeP1_none_1.normr;

%validMeP2_none_2 = ~isnan(MeP2_none_2);
%validdata_none_2 = ~isnan(data_none_2);
%validdataBoth_none_2 = validMeP2_none_2 & validdata_none_2;
%MeP2_none_valid_2 = MeP2_none_2(validdataBoth_none_2) ;
%data_valid_none_2 = data_none_2(validdataBoth_none_2);

%[p_MeP2_none_2,e_MeP2_none_2]=polyfit(MeP2_none_valid_2, data_valid_none_2,1);
%err_MeP2_none_2=e_MeP2_none_2.normr;

%validMeP1_ifn1_1 = ~isnan(MeP1_ifn1_1);
%validdata_ifn1_1 = ~isnan(data_ifn1_1);
%validdataBoth_ifn1_1 = validMeP1_ifn1_1 & validdata_ifn1_1;
%MeP1_ifn1_valid_1 = MeP1_ifn1_1(validdataBoth_ifn1_1) ;
%data_valid_ifn1_1 = data_ifn1_1(validdataBoth_ifn1_1);

%[p_MeP1_ifn1_1,e_MeP1_ifn1_1]=polyfit(MeP1_ifn1_valid_1, data_valid_ifn1_1,1);
%err_MeP1_ifn1_1=e_MeP1_ifn1_1.normr;

%validMeP2_ifn1_2 = ~isnan(MeP2_ifn1_2);
%validdata_ifn1_2 = ~isnan(data_ifn1_2);
%validdataBoth_ifn1_2 = validMeP2_ifn1_2 & validdata_ifn1_2;
%MeP2_ifn1_valid_2 = MeP2_ifn1_2(validdataBoth_ifn1_2) ;
%data_valid_ifn1_2 = data_ifn1_2(validdataBoth_ifn1_2);

%[p_MeP2_ifn1_2,e_MeP2_ifn1_2]=polyfit(MeP2_ifn1_valid_2, data_valid_ifn1_2,1);
%err_MeP2_ifn1_2=e_MeP2_ifn1_2.normr;


%err=err_MeP1_none_1+err_MeP2_none_2 + err_MeP1_ifn1_1 + err_MeP2_ifn1_2;

validMeP1_both_none_1 = ~isnan(MeP1_none_1);
validMeP1_both_ifn1_1 = ~isnan(MeP1_ifn1_1);
validdata_both_none_1 = ~isnan(data_none_1);
validdata_both_ifn1_1 = ~isnan(data_ifn1_1);
validdataAll_both_1 = validMeP1_both_none_1 &validMeP1_both_ifn1_1 & validdata_both_none_1 & validdata_both_ifn1_1 ;
%MeP1_both_valid_1 = MeP1_none_1(validdataAll_none_1) ;
%data_valid_none_1 = data_none_1(validdataBoth_none_1);

%min_MeP1_none_1 = min(MeP1_none_1(validdataAll_both_1));
%min_MeP1_ifn1_1 = min(MeP1_ifn1_1(validdataAll_both_1));
%max_MeP1_none_1 = max(MeP1_none_1(validdataAll_both_1));
%max_MeP1_ifn1_1 = max(MeP1_ifn1_1(validdataAll_both_1));

%norm_MeP1_none_1 = (MeP1_none_1(validdataAll_both_1) - min_MeP1_none_1) / ( max_MeP1_none_1 - min_MeP1_none_1); 
%your_original_data = minVal + norm_data.*(maxVal - minVal)
%norm_MeP1_ifn1_1 = (MeP1_ifn1_1(validdataAll_both_1) - min_MeP1_ifn1_1) / ( max_MeP1_ifn1_1 - min_MeP1_ifn1_1);

MeP1_both=[MeP1_none_1(validdataAll_both_1) , MeP1_ifn1_1(validdataAll_both_1) ];
%MeP1_both=[norm_MeP1_none_1 , norm_MeP1_ifn1_1 ];

%min_data_none_1 = min(data_none_1(validdataAll_both_1));
%%min_data_ifn1_1 = min(data_ifn1_1(validdataAll_both_1));
%max_data_none_1 = max(data_none_1(validdataAll_both_1));
%max_data_ifn1_1 = max(data_ifn1_1(validdataAll_both_1));

%norm_data_none_1 = (data_none_1(validdataAll_both_1) - min_data_none_1) / ( max_data_none_1 - min_data_none_1); 
%your_original_data = minVal + norm_data.*(maxVal - minVal)
%norm_data_ifn1_1 = (data_ifn1_1(validdataAll_both_1) - min_data_ifn1_1) / ( max_MeP1_ifn1_1 - min_MeP1_ifn1_1);

data_MeP1_both=[data_none_1(validdataAll_both_1), data_ifn1_1(validdataAll_both_1)];

%[p_MeP1_both, e_MeP1_both]=polyfit(MeP1_both, data_MeP1_both, 1);

err_MeP1_none_1 = sum((sf3*MeP1_none_1(validdataAll_both_1) - data_none_1(validdataAll_both_1)).^2);
err_MeP1_ifn1_1 = sum((sf3*MeP1_ifn1_1(validdataAll_both_1) - data_ifn1_1(validdataAll_both_1)).^2);


%validMeP2_both_none_2 = ~isnan(MeP2_none_2);
%validMeP2_both_ifn1_2 = ~isnan(MeP2_ifn1_2);
%validdata_both_none_2 = ~isnan(data_none_2);
%validdata_both_ifn1_2 = ~isnan(data_ifn1_2);
%validdataAll_both_2 = validMeP2_both_none_2 &validMeP2_both_ifn1_2 & validdata_both_none_2 & validdata_both_ifn1_2 ;

%min_MeP2_none_2 = min(MeP2_none_2(validdataAll_both_2));
%min_MeP2_ifn1_2 = min(MeP2_ifn1_2(validdataAll_both_2));
%max_MeP2_none_2 = max(MeP2_none_2(validdataAll_both_2));
%max_MeP2_ifn1_2 = max(MeP2_ifn1_2(validdataAll_both_2));

%norm_MeP2_none_2 = (MeP2_none_2(validdataAll_both_2) - min_MeP2_none_2) / ( max_MeP2_none_2 - min_MeP2_none_2); 
%your_original_data = minVal + norm_data.*(maxVal - minVal)
%norm_MeP2_ifn1_2 = (MeP2_ifn1_2(validdataAll_both_2) - min_MeP2_ifn1_2) / ( max_MeP2_ifn1_2 - min_MeP2_ifn1_2);

%MeP2_both=[MeP2_none_2(validdataAll_both_2) , MeP2_ifn1_2(validdataAll_both_2)];
%MeP2_both=[norm_MeP2_none_2 , norm_MeP2_ifn1_2];

%min_data_none_2 = min(data_none_2(validdataAll_both_2));
%min_data_ifn1_2 = min(data_ifn1_2(validdataAll_both_2));
%max_data_none_2 = max(data_none_2(validdataAll_both_2));
%max_data_ifn1_2 = max(data_ifn1_2(validdataAll_both_2));

%norm_data_none_2 = (data_none_2(validdataAll_both_2) - min_data_none_2) / ( max_data_none_2 - min_data_none_2); 
%your_original_data = minVal + norm_data.*(maxVal - minVal)
%norm_data_ifn1_2 = (data_ifn1_2(validdataAll_both_2) - min_data_ifn1_2) / ( max_data_ifn1_2 - min_data_ifn1_2);

%data_MeP2_both=[data_none_2(validdataAll_both_2), data_ifn1_2(validdataAll_both_2)];
%[p_MeP2_both, e_MeP2_both]=polyfit(MeP2_both, data_MeP2_both, 1);

%err_MeP2_none_2 = sum((sf4*MeP2_none_2(validdataAll_both_2) - data_none_2(validdataAll_both_2)).^2);
%err_MeP2_ifn1_2 = sum((sf4*MeP2_ifn1_2(validdataAll_both_2) - data_ifn1_2(validdataAll_both_2)).^2);

%err_MeP1_both = e_MeP1_both.normr;
%err_MeP2_both = e_MeP2_both.normr;

err=err_MeP1_none_1 + err_MeP1_ifn1_1;%+ err_MeP2_none_2 + err_MeP2_ifn1_2; %
sf
p_MeP1_none_1=[sf3 0];%;%p_MeP1_both;%
p_MeP2_none_2=[sf3 0];%p_MeP2_both;%
p_MeP1_ifn1_1=[sf3 0];%p_MeP1_both;%
p_MeP2_ifn1_2=[sf3 0];%p_MeP2_both;%

%p_MeP1_both
%p_MeP2_both
cmap=hsv(8);
figure(8);
semilogx(gtarget_none_1(:,1),p_MeP1_none_1(1)*MeP1_none_1(:,1)+p_MeP1_none_1(2),gtarget_none_1(:,1),data_none_1(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_none_1(:,2),p_MeP1_none_1(1)*MeP1_none_1(:,2)+p_MeP1_none_1(2),gtarget_none_1(:,2),data_none_1(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,3),p_MeP1_none_1(1)*MeP1_none_1(:,3)+p_MeP1_none_1(2),gtarget_none_1(:,3),data_none_1(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,4),p_MeP1_none_1(1)*MeP1_none_1(:,4)+p_MeP1_none_1(2),gtarget_none_1(:,4),data_none_1(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,5),p_MeP1_none_1(1)*MeP1_none_1(:,5)+p_MeP1_none_1(2),gtarget_none_1(:,5),data_none_1(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,6),p_MeP1_none_1(1)*MeP1_none_1(:,6)+p_MeP1_none_1(2),gtarget_none_1(:,6),data_none_1(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,7),p_MeP1_none_1(1)*MeP1_none_1(:,7)+p_MeP1_none_1(2),gtarget_none_1(:,7),data_none_1(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,8),p_MeP1_none_1(1)*MeP1_none_1(:,8)+p_MeP1_none_1(2),gtarget_none_1(:,8),data_none_1(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
figure(9);
%subplot(1,2,2)
semilogx(gtarget_none_2(:,1),p_MeP2_none_2(1)*MeP2_none_2(:,1)+p_MeP2_none_2(2),gtarget_none_2(:,1),data_none_2(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_none_2(:,2),p_MeP2_none_2(1)*MeP2_none_2(:,2)+p_MeP2_none_2(2),gtarget_none_2(:,2),data_none_2(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,3),p_MeP2_none_2(1)*MeP2_none_2(:,3)+p_MeP2_none_2(2),gtarget_none_2(:,3),data_none_2(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,4),p_MeP2_none_2(1)*MeP2_none_2(:,4)+p_MeP2_none_2(2),gtarget_none_2(:,4),data_none_2(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,5),p_MeP2_none_2(1)*MeP2_none_2(:,5)+p_MeP2_none_2(2),gtarget_none_2(:,5),data_none_2(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,6),p_MeP2_none_2(1)*MeP2_none_2(:,6)+p_MeP2_none_2(2),gtarget_none_2(:,6),data_none_2(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,7),p_MeP2_none_2(1)*MeP2_none_2(:,7)+p_MeP2_none_2(2),gtarget_none_2(:,7),data_none_2(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_2(:,8),p_MeP2_none_2(1)*MeP2_none_2(:,8)+p_MeP2_none_2(2),gtarget_none_2(:,8),data_none_2(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off

figure(10);
semilogx(gtarget_ifn1_1(:,1),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,1)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,1),data_ifn1_1(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_ifn1_1(:,2),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,2)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,2),data_ifn1_1(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,3),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,3)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,3),data_ifn1_1(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,4),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,4)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,4),data_ifn1_1(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,5),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,5)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,5),data_ifn1_1(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,6),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,6)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,6),data_ifn1_1(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,7),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,7)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,7),data_ifn1_1(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,8),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,8)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,8),data_ifn1_1(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
figure(11);
semilogx(gtarget_ifn1_2(:,1),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,1)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,1),data_ifn1_2(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_ifn1_2(:,2),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,2)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,2),data_ifn1_2(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,3),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,3)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,3),data_ifn1_2(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,4),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,4)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,4),data_ifn1_2(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,5),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,5)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,5),data_ifn1_2(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,6),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,6)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,6),data_ifn1_2(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,7),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,7)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,7),data_ifn1_2(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_2(:,8),p_MeP2_ifn1_2(1)*MeP2_ifn1_2(:,8)+p_MeP2_ifn1_2(2),gtarget_ifn1_2(:,8),data_ifn1_2(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off
return

