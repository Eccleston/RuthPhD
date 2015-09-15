%% Process CliLBS output

clear
close all

fname = 'check6_gMfixed';
trial = 1;

doSimulation = 1;
doParameters = 1;

f.saving = 1;
f.dataDir = '../CliLBS/';
outname = sprintf('%s_tr%d',fname,trial);
f.simDir = sprintf('../CliLBS/%s_inference/Trial%d/',fname,trial);
if doSimulation
  f.locs = [8 8 8 8 8 8 8 8];
  f.title = 'No IFN';
  f.outname = [outname '_none'];
  fh1 = ProcessOutput('SSL_surf_none_data','mysweep1',f);
  f.locs = [8 8 7 8 7 8 8 8];
  f.title = '1 ng IFN';
  f.outname = [outname '_ifn'];
  fh2 = ProcessOutput('SSL_surf_IFN1_data','mysweep2',f);
  fh2.Position = [100 100 450 300];
end
if doParameters
  f.outname = outname;
  post = ProcessParameters(f,{});
end

return