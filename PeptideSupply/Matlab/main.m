%% Main file for modelling pMHC abundance

clear all
warning off

%% Load data

datadir = '../Data/';

dat1 = dlmread([datadir '150401_ASN_SSL.txt'],'\t',1,1);
locs = 3:66;
in1_none = reshape(dat1(locs,1),8,8);
target1_none = reshape(dat1(locs,2),8,8);
in1_ifn1 = reshape(dat1(locs,3),8,8);
target1_ifn1 = reshape(dat1(locs,4),8,8);

dat2 = dlmread([datadir '150422_ASN_SSL.txt'],'\t',1,1);
locs = 3:66;
in2_none = reshape(dat2(locs,1),8,8);
target2_none = reshape(dat2(locs,2),8,8);
in2_ifn1 = reshape(dat2(locs,3),8,8);
target2_ifn1 = reshape(dat2(locs,4),8,8);
in2_ifn2 = reshape(dat2(locs,5),8,8);
target2_ifn2 = reshape(dat2(locs,6),8,8);

%% Simulate the model (test)
u1 = 1.6e-5;
u2 = 3.4e-5;
gs = 10.^(0:0.5:4);
Ng = length(gs);
tic
for i1 = 1:Ng
  g1 = gs(i1);
  for i2 = 1:Ng
    g2 = gs(i2);
    [MeP1(i1,i2),MeP2(i1,i2)] = simulateMHC(g1,g2,u1,u2);
  end
end
toc

figure(1)
subplot(1,2,1)
semilogx(gs,MeP1)
subplot(1,2,2)
semilogx(gs,MeP2)

%% Fit the scale factors

% Set g1 to be the measurement of the target peptide (cytoplasm)
% Set g2 to be the measurement of the competitor peptide (cytoplasm)
% Set data to be the cell surface measure of the target peptide

%sfOpt = fminsearch(@(sf) sum((sf(2)*simulateMHC(sf(1)*g1,sf(1)*g2,u1,u2) - data).^2),[1;1]);

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