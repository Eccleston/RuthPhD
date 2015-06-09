%% Create a plot of b5703_short results

Ns = 10;
workingdir = 'b5703_short_simulation/';

for i= 1:Ns
  data = dlmread([workingdir 'stochastic' num2str(i) '.tsv'],'\t',1,0);
  time = data(:,1);
  SP = data(:,5:4:17);
  
  subplot(3,4,i)
  plot(time,SP)
  legend('SP4','SP3','SP2','SP1')
end

return