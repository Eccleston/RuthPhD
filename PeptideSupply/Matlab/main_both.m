%% Main file for modelling pMHC abundance

clear all
warning off

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

%dat1 = dlmread([datadir 'AllData_150401_ASN_SSL.txt'],' ',1,1);
%locs = 3:66;
%in1_none = reshape(dat1(locs,1),8,8);
%target1_none = reshape(dat1(locs,2),8,8);
%in1_ifn1 = reshape(dat1(locs,3),8,8);
%target1_ifn1 = reshape(dat1(locs,4),8,8);
%competitor1_none = reshape(dat1(locs,5),8,8);
%competitor1_ifn1 = reshape(dat1(locs,6),8,8);


%dat2 = dlmread([datadir 'AllData_150422_ASN_SSL.txt'],' ',1,1);
%locs = 3:66;
%in2_none = reshape(dat2(locs,1),8,8);
%target2_none = reshape(dat2(locs,2),8,8);
%in2_ifn1 = reshape(dat2(locs,3),8,8);
%target2_ifn1 = reshape(dat2(locs,4),8,8);
%in2_ifn2 = reshape(dat2(locs,5),8,8);
%target2_ifn2 = reshape(dat2(locs,6),8,8);
%competitor2_none = reshape(dat2(locs,7),8,8);
%competitor2_ifn1 = reshape(dat2(locs,8),8,8);
%competitor2_ifn2 = reshape(dat2(locs,9),8,8);

%dat3 = dlmread([datadir 'AllData_150520_ASN_SSL.txt'],',',1,1);
%locs = 3:66;
%SSL_cyt_none = reshape(dat3(locs,1),8,8);
%SSL_surf_none = reshape(dat3(locs,2),8,8);
%SSL_cyt_ifn1 = reshape(dat3(locs,3),8,8);
%SSL_surf_ifn1 = reshape(dat3(locs,4),8,8);
%SSL_cyt_ifn2 = reshape(dat3(locs,5),8,8);
%SSL_surf_ifn2 = reshape(dat3(locs,6),8,8);
%ASN_cyt_none = reshape(dat3(locs,7),8,8);
%ASN_cyt_ifn1 = reshape(dat3(locs,8),8,8);
%ASN_cyt_ifn2 = reshape(dat3(locs,9),8,8);%

%dat4 = dlmread([datadir 'AllData_150520_ASN_SSL_reverse.txt'],',',1,1);
%locs = 3:66;
%SSL_cyt_none_2 = reshape(dat4(locs,1),8,8);
%ASN_surf_none = reshape(dat4(locs,2),8,8);
%SSL_cyt_ifn1_2 = reshape(dat4(locs,3),8,8);
%ASN_surf_ifn1 = reshape(dat4(locs,4),8,8);
%SSL_cyt_ifn2_2 = reshape(dat4(locs,5),8,8);
%ASN_surf_ifn2 = reshape(dat4(locs,6),8,8);
%ASN_cyt_none_2 = reshape(dat4(locs,7),8,8);
%ASN_cyt_ifn1_2 = reshape(dat4(locs,8),8,8);
%ASN_cyt_ifn2_2 = reshape(dat4(locs,9),8,8);

%% Simulate the model (test)

% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);
% Set g1 to be the measurement of the target peptide (cytoplasm)
gtarget_none=SSL_cyt_none;%in2_none;
% Set g2 to be the measurement of the competitor peptide (cytoplasm)
gcomp_none=ASN_cyt_none;%competitor2_none;

Ng=8;
gM=150.5;
for i1 = 1:Ng
    for i2 = 1:Ng
    
  g1 = gtarget_none(i1,i2); % 
  g2 = gcomp_none(i1,i2);
    [MeP11(i1,i2),MeP22(i1,i2)] = simulateMHC(g1,g2,u1,u2,gM);
  end
end

msize = 10;
figure(1)
subplot(1,2,1)
semilogx(gtarget_none,MeP11);legend({'gcomp1','gcomp2','gcomp3','gcomp4','gcomp5','gcomp6','gcomp7','gcomp8'})
%axis([xlims ylims1])
subplot(1,2,2)
semilogx(gtarget_none,MeP22)
%axis([xlims ylims2])

%% Fit the scale factors
% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);
% Set g1_none to be the measurement of the target peptide (cytoplasm) 
gtarget_none=SSL_cyt_none;%
% Set g2 to be the measurement of the competitor peptide (cytoplasm) 
gcomp_none=ASN_cyt_none;
% Set data to be the cell surface measure of the target peptide
data_none = SSL_surf_none;
% Set g1 to be the measurement of the target peptide (cytoplasm) *sf1
gtarget_ifn=SSL_cyt_ifn1;%
% Set g2 to be the measurement of the competitor peptide (cytoplasm) *sf2
gcomp_ifn=ASN_cyt_ifn1;
% Set data to be the cell surface measure of the target peptide
data_ifn = SSL_surf_ifn1;
sf=[0.1323,0.78,0.00039,1537];%[1,1,105.5,10000.5];

[sf l_s]=fminsearch(@(sf)least_squares_both(sf,gtarget_none,gcomp_none,gtarget_ifn,gcomp_ifn,data_none,data_ifn),sf);


%% Plot
xlims = [10 1e6];
ylims1 = [0 14000];
ylims2 = [0 100000];
msize = 10;

figure(2);
subplot(2,3,1)
set(gca,'NextPlot','ReplaceChildren','ColorOrder',hsv(10))
semilogx(in1_none,target1_none,'.-','MarkerSize',msize)
axis([xlims ylims1])
box off
subplot(2,3,2)
semilogx(in1_ifn1,target1_ifn1,'.-','MarkerSize',msize)
axis([xlims ylims2])
box off
%subplot(2,3,3)
%semilogx(gtarget,*MeP11)
%axis([xlims ylims1])
%box off
subplot(2,3,4)
semilogx(in2_none,target2_none,'.-','MarkerSize',msize)
axis([xlims ylims1])
box off
subplot(2,3,5)
semilogx(in2_ifn1,target2_ifn1,'.-','MarkerSize',msize)
axis([xlims ylims2])
box off
subplot(2,3,6)
semilogx(in2_ifn2,target2_ifn2,'.-','MarkerSize',msize)
axis([xlims ylims2])
box off



return