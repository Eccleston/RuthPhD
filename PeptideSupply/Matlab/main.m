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

%% Simulate the model

%% Plot
xlims = [10 1e6];
ylims1 = [0 14000];
ylims2 = [0 100000];

subplot(2,3,1)
%set(gca,'NextPlot','ReplaceChildren','ColorOrder',jet(7))
semilogx(in1_none,target1_none,'o-')
axis([xlims ylims1])

box off
subplot(2,3,2)
semilogx(in1_ifn1,target1_ifn1,'o-')
axis([xlims ylims2])
box off

subplot(2,3,4)
semilogx(in2_none,target2_none,'o-')
axis([xlims ylims1])
box off
subplot(2,3,5)
semilogx(in2_ifn1,target2_ifn1,'o-')
axis([xlims ylims2])
box off
subplot(2,3,6)
semilogx(in2_ifn2,target2_ifn2,'o-')
axis([xlims ylims2])
box off

return