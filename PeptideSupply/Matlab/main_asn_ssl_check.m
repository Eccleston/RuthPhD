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
%% Simulate the model (test)

%u1 = 1.6e-5;
%u2 = 3.4e-5;
%gs = 10.^(0:0.5:4);
%Ng = length(gs);

% Set u1 to be the measurement of the target peptide off-rate
u1 = log(2)/(728.5746 * 60);
% Set u2 to be the measurement of the competitor peptide off-rate
u2 = log(2)/(337.9968 * 60);
% Set g1 to be the measurement of the target peptide (cytoplasm)
gtarget=SSL_cyt_none_2;
% Set g2 to be the measurement of the competitor peptide (cytoplasm)
gcomp=ASN_cyt_none_2;
%tfinal=1*3600;
%tf=24*3600;
Ng=8;

for i1 = 1:Ng
    for i2 = 1:Ng
   
  g1 = gtarget(i1,i2); %
 
   g2 = gcomp(i1,i2);
    [MeP11(i1,i2),MeP22(i1,i2)] = simulateMHC(g1,g2,u1,u2,150.5);
  end
end

msize = 10;
figure(1)
subplot(1,2,1)
semilogx(gtarget,MeP11);legend({'gcomp1','gcomp2','gcomp3','gcomp4','gcomp5','gcomp6','gcomp7','gcomp8'})
%axis([xlims ylims1])
subplot(1,2,2)
semilogx(gtarget,MeP22)
%axis([xlims ylims2])

%% Fit the scale factors
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
gtarget_ifn1_1=SSL_cyt_ifn1;
% set g2_2 to be the second measurement of the competitor peptide
% (cytoplasm)
gcomp_ifn1_1=ASN_cyt_ifn1;
% Set data_2 to be the cell surface measure of the competitor peptide
data_ifn1_1=SSL_surf_ifn1;
% Set g1 to be the measurement of the target peptide (cytoplasm) 
gtarget_ifn1_2=SSL_cyt_ifn1_2;%
% Set g2 to be the measurement of the competitor peptide (cytoplasm) 
gcomp_ifn1_2=ASN_cyt_ifn1_2;
% Set data to be the cell surface measure of the target peptide
data_ifn1_2 = ASN_surf_ifn1;


sf=[1,1,1,1];

[sf l_s]=fminsearch(@(sf)least_squares_asn_ssl_ifn_check(sf,gtarget_none_1,gcomp_none_1,data_none_1,gtarget_none_2,gcomp_none_2, data_none_2,gtarget_ifn1_1,gcomp_ifn1_1,data_ifn1_1,gtarget_ifn1_2,gcomp_ifn1_2, data_ifn1_2),sf);

for i1 = 1:Ng
    for i2 = 1:Ng
    %g1=gs(i1);
  g1 = gtarget(i1,i2); % 
  %    for i2 = 1:Ng
%g2=gs(i2);
   g2 = gcomp(i1,i2);
    [MeP11(i1,i2),MeP22(i1,i2),M(i1,i2),T(i1,i2),MT(i1,i2)] = simulateMHC(sf(1)*g1,sf(2)*g2,u1,u2,150.5);
  end
end


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