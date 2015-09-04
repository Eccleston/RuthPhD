
sf1=0.1323;%0.109191165981564;%0.00504133178459042;%0.003411838;

sf2=0.7804;%0.555429420251379;%0.969957735248499;%0.237878986;

sf3=3.9468e-4;%0.00100010131001115;%0.00131976368681847;%0.004997214;

gM=150.5;

upfactor=1.5371e3;%390.7520171688;%0;%999.228892;

%% Load data

datadir = '../Data/';

dat1 = dlmread([datadir 'AllData_150520_ASN_SSL.txt'],',',1,1);
locs = 3:66;
SSL_cyt_none = reshape(dat1(locs,1),8,8);
SSL_surf_none = reshape(dat1(locs,2),8,8);
SSL_cyt_ifn1 = reshape(dat1(locs,3),8,8);
SSL_surf_ifn1 = reshape(dat1(locs,4),8,8);
SSL_cyt_ifn2 = reshape(dat1(locs,5),8,8);
SSL_surf_ifn2 = reshape(dat1(locs,6),8,8);
ASN_cyt_none = reshape(dat1(locs,7),8,8);
ASN_cyt_ifn1 = reshape(dat1(locs,8),8,8);
ASN_cyt_ifn2 = reshape(dat1(locs,9),8,8);

% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);

% Set gtarget_none to be the measurement of the target peptide
% (cytoplasm)without IFN
gtarget_none_1=SSL_cyt_none;%
% Set gcomp_none to be the measurement of the competitor peptide
% (cytoplasm) without IFN
gcomp_none_1=ASN_cyt_none;
% Set data_none to be the cell surface measure of the target peptide
% without IFN
data_none_1 = SSL_surf_none;

% Set g1_2 to be the second measurement of the target peptide (cytoplasm)
gtarget_ifn1_1=SSL_cyt_ifn1(1:7,1:7);
% set g2_2 to be the second measurement of the competitor peptide
% (cytoplasm)
gcomp_ifn1_1=ASN_cyt_ifn1(1:7,1:7);
% Set data_2 to be the cell surface measure of the competitor peptide
data_ifn1_1=SSL_surf_ifn1(1:7,1:7);

Ng1=8; Ng2=8;

for i1=1:Ng1
    for i2=1:Ng1
                
        [MeP1_none_1(i1, i2), MeP2_none_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_none_1(i1, i2), sf2*gcomp_none_1(i1, i2), u1, u2, gM, 0, 0);
        %[MeP1_none_2(i1, i2), MeP2_none_2(i1, i2)]=simulateMHC_asn_ssl(sf1*gtarget_none_2(i1, i2), sf2*gcomp_none_2(i1, i2), u1, u2, gM1);
        %[MeP1_ifn1_1(i1, i2), MeP2_ifn1_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_1(i1, i2), sf2*gcomp_ifn1_1(i1, i2), u1, u2, gM, 1, upfactor);
        %[MeP1_ifn1_2(i1, i2), MeP2_ifn1_2(i1, i2)]=simulateMHC_asn_ssl(sf1*gtarget_ifn1_2(i1, i2), sf2*gcomp_ifn1_2(i1, i2), u1, u2, upreg*gM1);
    end
end
for i1=1:7
    for i2=1:7
                
        %[MeP1_none_1(i1, i2), MeP2_none_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_none_1(i1, i2), sf2*gcomp_none_1(i1, i2), u1, u2, gM, 0, upfactor);
        %[MeP1_none_2(i1, i2), MeP2_none_2(i1, i2)]=simulateMHC_asn_ssl(sf1*gtarget_none_2(i1, i2), sf2*gcomp_none_2(i1, i2), u1, u2, gM1);
        [MeP1_ifn1_1(i1, i2), MeP2_ifn1_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_1(i1, i2), sf2*gcomp_ifn1_1(i1, i2), u1, u2, gM, 1, upfactor);
        %[MeP1_ifn1_2(i1, i2), MeP2_ifn1_2(i1, i2)]=simulateMHC_asn_ssl(sf1*gtarget_ifn1_2(i1, i2), sf2*gcomp_ifn1_2(i1, i2), u1, u2, upreg*gM1);
   end
end
%validMeP1_both_none_1 = ~isnan(MeP1_none_1);
%validMeP1_both_ifn1_1 = ~isnan(MeP1_ifn1_1);
%validdata_both_none_1 = ~isnan(data_none_1);
%validdata_both_ifn1_1 = ~isnan(data_ifn1_1);
%validdataAll_both_1 = validMeP1_both_none_1 &validMeP1_both_ifn1_1 & validdata_both_none_1 & validdata_both_ifn1_1 ;
p_MeP1_none_1=[sf3 0];%;%p_MeP1_both;
%p_MeP2_none_2=[sf4 0];%p_MeP2_both;
p_MeP1_ifn1_1=[sf3 0];%p_MeP1_both;
%p_MeP2_ifn1_2=[sf4 0];%p_MeP2_both;

%p_MeP1_both
%p_MeP2_both
cmap=hsv(8);
figure(10);
semilogx(gtarget_none_1(:,1),sf3*MeP1_none_1(:,1),gtarget_none_1(:,1),data_none_1(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_none_1(:,2),sf3*MeP1_none_1(:,2),gtarget_none_1(:,2),data_none_1(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,3),sf3*MeP1_none_1(:,3),gtarget_none_1(:,3),data_none_1(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,4),sf3*MeP1_none_1(:,4),gtarget_none_1(:,4),data_none_1(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,5),sf3*MeP1_none_1(:,5),gtarget_none_1(:,5),data_none_1(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,6),sf3*MeP1_none_1(:,6),gtarget_none_1(:,6),data_none_1(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,7),sf3*MeP1_none_1(:,7),gtarget_none_1(:,7),data_none_1(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_none_1(:,8),sf3*MeP1_none_1(:,8),gtarget_none_1(:,8),data_none_1(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off

figure(11);
semilogx(gtarget_ifn1_1(:,1),sf3*MeP1_ifn1_1(:,1),gtarget_ifn1_1(:,1),data_ifn1_1(:,1),'o','Color',cmap(1,:),'LineWidth',2,'MarkerSize',10);hold on;
semilogx(gtarget_ifn1_1(:,2),sf3*MeP1_ifn1_1(:,2),gtarget_ifn1_1(:,2),data_ifn1_1(:,2),'o','Color',cmap(2,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,3),sf3*MeP1_ifn1_1(:,3),gtarget_ifn1_1(:,3),data_ifn1_1(:,3),'o','Color',cmap(3,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,4),sf3*MeP1_ifn1_1(:,4),gtarget_ifn1_1(:,4),data_ifn1_1(:,4),'o','Color',cmap(4,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,5),sf3*MeP1_ifn1_1(:,5),gtarget_ifn1_1(:,5),data_ifn1_1(:,5),'o','Color',cmap(5,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,6),sf3*MeP1_ifn1_1(:,6),gtarget_ifn1_1(:,6),data_ifn1_1(:,6),'o','Color',cmap(6,:),'LineWidth',2,'MarkerSize',10);
semilogx(gtarget_ifn1_1(:,7),sf3*MeP1_ifn1_1(:,7),gtarget_ifn1_1(:,7),data_ifn1_1(:,7),'o','Color',cmap(7,:),'LineWidth',2,'MarkerSize',10);
%semilogx(gtarget_ifn1_1(:,8),p_MeP1_ifn1_1(1)*MeP1_ifn1_1(:,8)+p_MeP1_ifn1_1(2),gtarget_ifn1_1(:,8),data_ifn1_1(:,8),'o','Color',cmap(8,:),'LineWidth',2,'MarkerSize',10);
xlabel('mCherry MFI','FontSize',16);ylabel('1C3 + GAM-AF647','FontSize',16);
legend({'P2 P9',' ','P10 P17',' ','P18 P25',' ','P26 P33',' ','P34 P40',' ','P42 P39',' ','P50 P57',' ','P58 P65',' '})
%legend({'V-ASN P2',' ','V-ASN P10',' ','V-ASN P18',' ','V-ASN P26',' ','V-ASN P34 ',' ','V-ASN P42',' ','V-ASN P50',' ','V-ASN P58',' '})
set(gca, 'FontSize',16);
hold off