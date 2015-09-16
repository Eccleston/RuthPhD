%% Process CliLBS output

clear
close all


%fname = 'check6_gMfixed'; trial = 1;
fname = 'check11_nd'; trial = 1;

outname = sprintf('%s_tr%d',fname,trial);
loadCases;

%cases = [SSL_none SSL_ifn];
cases = [SSL_none SSL_ifn ASN_none ASN_ifn];

doSimulation = 1;
doParameters = 1;

f.saving = 1;
f.dataDir = '../CliLBS/';
f.simDir = sprintf('../CliLBS/%s_inference/Trial%d/',fname,trial);
if doSimulation
  for i = 1:length(cases)
    f.locs = cases(i).locs;
    f.title = cases(i).title;
    f.outname = cases(i).outname;
    fh(i) = ProcessOutput(cases(i).datafile,cases(i).sweep,f);
  end
end
if doParameters
  f.outname = outname;
  post = ProcessParameters(f,{});
end

return