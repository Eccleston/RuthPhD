%% Main file for modelling pMHC abundance

clear all
warning off

%% Load data

datadir = '../Data/';

dat1 = dlmread([datadir 'AllData_150401_ASN_SSL.txt'],' ',1,1);
locs = 3:66;
in1_none = reshape(dat1(locs,1),8,8);
target1_none = reshape(dat1(locs,2),8,8);
in1_ifn1 = reshape(dat1(locs,3),8,8);
target1_ifn1 = reshape(dat1(locs,4),8,8);
competitor1_none = reshape(dat1(locs,5),8,8);
competitor1_ifn1 = reshape(dat1(locs,6),8,8);


dat2 = dlmread([datadir 'AllData_150422_ASN_SSL.txt'],' ',1,1);
locs = 3:66;
in2_none = reshape(dat2(locs,1),8,8);
target2_none = reshape(dat2(locs,2),8,8);
in2_ifn1 = reshape(dat2(locs,3),8,8);
target2_ifn1 = reshape(dat2(locs,4),8,8);
in2_ifn2 = reshape(dat2(locs,5),8,8);
target2_ifn2 = reshape(dat2(locs,6),8,8);
competitor2_none = reshape(dat2(locs,7),8,8);
competitor2_ifn1 = reshape(dat2(locs,8),8,8);
competitor2_ifn2 = reshape(dat2(locs,9),8,8);
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
gtarget=in2_none;
% Set g2 to be the measurement of the competitor peptide (cytoplasm)
gcomp=competitor2_none;
%tfinal=1*3600;
%tf=24*3600;
Ng=8;

for i1 = 1:Ng
    for i2 = 1:Ng
    %g1=gs(i1);
  g1 = gtarget(i1,i2); % 
  %    for i2 = 1:Ng
%g2=gs(i2);
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
% Set g1 to be the measurement of the target peptide (cytoplasm) *sf1

gtarget=in2_none;%
% Set g2 to be the measurement of the competitor peptide (cytoplasm) *sf2

gcomp=competitor2_none;
% Set data to be the cell surface measure of the target peptide
data = target2_none;
%fun=@(sf)sum((sf(2)*simulateMHC(sf(1)*g1,sf(1)*g2,u1,u2)-data).^2)
%sfOpt = fminsearch(least_square_fun,[1;1]);
sf=[10,10];
%Estimate paramters and variances
%LB=[0, 0, 0];%60*0.00019254];%Lower bounds of every parameter
%%UB=[1e5,1e5,1e5];%500*60*0.00019254];%Upper bounds of every parameter
%optim_options=optimset('Display', 'iter',...
%    'TolFun',1e-8,...
%    'TolX', 1e-8,...
%    'MaxIter',1e9,...%Maximum number of iterations allowed.
%    'MaxFunEvals',1E4,...%Maximum number of function evaluations allowed
%    'Algorithm','levenberg-marquardt'); 

%[k,resnorm, residual, exitflag, OUTPUT, LAMBDA, Jacobian]=lsqnonlin(@least_squares_fun, sf0,[] ,[],optim_options,gtarget,gcomp,data,target1_none);
%%disp(' ')
%sf
[sf l_s]=fminsearch(@(sf)least_squares_poly(sf,gtarget,gcomp,data),sf);
%upreg=sf;
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