%% Process the inference from LBS
function fh = ProcessOutput(dataFile,simFile,f)
%dataFile = 'SSL_surf_none_data'; simFile = 'mysweep1';

[hd,datI] = importTSV([f.dataDir dataFile '.txt']);
simI = dlmread([f.simDir simFile '.tsv'],'\t'); simI = simI(2:2:end);
k = 0;
for i = 1:length(f.locs)
  for j = 1:f.locs(i)
    str = strsplit(hd{k+j+1},'_');
    supp1{i}(j) = str2double(str{2});
    supp2{i}(j) = str2double(str{3});
  end
  data{i} = datI(k+2:k+f.locs(i)+1);
  sim{i} = simI(k+1:k+f.locs(i));
  k = k + f.locs(i);
end
N = size(sim,2);
A = hsv(sum(N));

%% Compute R^2
%plot(log10(datI(2:end)),log10(simI),'x')
R2 = corrcoef(log10(datI(2:end)),log10(simI));
r = R2(1,2);

%% Compare model MLE versus data
fh = figure;
fh.Position = [100 500 450 300];
clf('reset')
hold on
for i = 1:8
  plot(supp1{i},data{i}/1000,'.','Color',A(i,:),'MarkerSize',15)
  plot(supp1{i},sim{i}/1000,'-','Color',A(i,:),'LineWidth',1.5)
end
set(gca,'Xscale','log')
hold off
xlabel('Input')
ylabel('Output (x 1,000)')
ylims = get(gca,'Ylim');
text(30,[0.1 0.9]*ylims',['\rho = ' num2str(r)])
title(f.title)

if f.saving
  figuresave(fh,[f.outname '_comparison'],300)
end

return