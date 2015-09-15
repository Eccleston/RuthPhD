
sf1=0.260196855793763;%0.265777100059837;%0.278388656956398;%0.260196855793763;%0.260943132268352;%0.382705104415568;%0.261653835936026;%0.0456116308509272;%0.00703887616435119;%0.194966;%0.109191165981564;%0.00504133178459042;%0.003411838;0.1323;%

sf2=1.04884576179473;%1.08752041751078;%0.0108695989165427;%1.04884576179473;%1.07378016444957;%5.20445828777391;%1.06351286346275;%0.329139724113906;%0.305983131587646;%0.5897358;%0.555429420251379;%0.969957735248499;%0.237878986;0.7804;%

sf3=0.0123838968437717;%0.0310572954053281;%1.10918299124869;%0.0123838968437717;%0.0127137308079439;%0.00567728148348992;%4.96872038305114;%0.00268859525507661;%0.0387020876481363;%0.00147364;%0.00100010131001115;%0.00131976368681847;%0.004997214;3.9468e-4;%

gM=10.751507530203;%5.75018046510033;%0.403712349349218;%10.751507530203;%10.4483443758914;%150.5;

upfactor=90.7471924699804;%38.5928862505608;%2.46768150345268;%90.7471924699804;%88.5524679633756;%0.0953102860792972;%0.00162947143327367;%872.53099829729;%0;%390.7520171688;%0;%999.228892;1.5371e3;%

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

dat2 = dlmread([datadir 'AllData_150520_ASN_SSL_reverse.txt'],',',1,1);
locs = 3:66;
SSL_cyt_none_2 = reshape(dat2(locs,1),8,8);
ASN_surf_none = reshape(dat2(locs,2),8,8);
SSL_cyt_ifn1_2 = reshape(dat2(locs,3),8,8);
ASN_surf_ifn1 = reshape(dat2(locs,4),8,8);
SSL_cyt_ifn2_2 = reshape(dat2(locs,5),8,8);
ASN_surf_ifn2 = reshape(dat2(locs,6),8,8);
ASN_cyt_none_2 = reshape(dat2(locs,7),8,8);
ASN_cyt_ifn1_2 = reshape(dat2(locs,8),8,8);
ASN_cyt_ifn2_2 = reshape(dat2(locs,9),8,8);

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

% Set gtarget_none to be the measurement of the target peptide
% (cytoplasm)without IFN
gtarget_none_2=SSL_cyt_none_2;%
% Set gcomp_none to be the measurement of the competitor peptide
% (cytoplasm) without IFN
gcomp_none_2=ASN_cyt_none_2;
% Set data_none to be the cell surface measure of the ASN peptide
% without IFN
data_none_2 = ASN_surf_none;
% Set g1_2 to be the second measurement of the target peptide (cytoplasm)
gtarget_ifn1_1=SSL_cyt_ifn1(1:7,:);
% set g2_2 to be the second measurement of the competitor peptide
% (cytoplasm)
gcomp_ifn1_1=ASN_cyt_ifn1(1:7,:);
% Set data_2 to be the cell surface measure of the competitor peptide
data_ifn1_1=SSL_surf_ifn1(1:7,:);
% Set g1 to be the measurement of the target peptide (cytoplasm) 
gtarget_ifn1_2=SSL_cyt_ifn1_2;%
% Set g2 to be the measurement of the competitor peptide (cytoplasm) 
gcomp_ifn1_2=ASN_cyt_ifn1_2;
% Set data to be the cell surface measure of the target peptide
data_ifn1_2 = ASN_surf_ifn1;


Ng1=8; Ng2=8;

for i1=1:Ng1
    for i2=1:Ng1
                
        [MeP1_none_1(i1, i2), MeP2_none_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_none_1(i1, i2), sf2*gcomp_none_1(i1, i2), u1, u2, gM, 0, 0);
        [MeP1_none_2(i1, i2), MeP2_none_2(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_none_2(i1, i2), sf2*gcomp_none_2(i1, i2), u1, u2, gM, 0, 0);
        %[MeP1_ifn1_1(i1, i2), MeP2_ifn1_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_1(i1, i2), sf2*gcomp_ifn1_1(i1, i2), u1, u2, gM, 1, upfactor);
        %[MeP1_ifn1_2(i1, i2), MeP2_ifn1_2(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_2(i1, i2), sf2*gcomp_ifn1_2(i1, i2), u1, u2, gM, 1, upfactor);
    end
end
for i1=1:7
    for i2=1:8
                
        %[MeP1_none_1(i1, i2), MeP2_none_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_none_1(i1, i2), sf2*gcomp_none_1(i1, i2), u1, u2, gM, 0, upfactor);
        %[MeP1_none_2(i1, i2), MeP2_none_2(i1, i2)]=simulateMHC_asn_ssl(sf1*gtarget_none_2(i1, i2), sf2*gcomp_none_2(i1, i2), u1, u2, gM1);
        [MeP1_ifn1_1(i1, i2), MeP2_ifn1_1(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_1(i1, i2), sf2*gcomp_ifn1_1(i1, i2), u1, u2, gM, 1, upfactor);
        [MeP1_ifn1_2(i1, i2), MeP2_ifn1_2(i1, i2)]=simulateMHC_mcmc(sf1*gtarget_ifn1_2(i1, i2), sf2*gcomp_ifn1_2(i1, i2), u1, u2, gM, 1, upfactor);
  end
end
%validMeP1_both_none_1 = ~isnan(MeP1_none_1);
%validMeP1_both_ifn1_1 = ~isnan(MeP1_ifn1_1);
%validdata_both_none_1 = ~isnan(data_none_1);
%validdata_both_ifn1_1 = ~isnan(data_ifn1_1);
%validdataAll_both_1 = validMeP1_both_none_1 &validMeP1_both_ifn1_1 & validdata_both_none_1 & validdata_both_ifn1_1 ;
p_MeP1_none_1=[sf3 0];%;%p_MeP1_both;
p_MeP2_none_2=[sf3 0];%p_MeP2_both;
p_MeP1_ifn1_1=[sf3 0];%p_MeP1_both;
p_MeP2_ifn1_2=[sf3 0];%p_MeP2_both;

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

figure(12);
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
figure(13);
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