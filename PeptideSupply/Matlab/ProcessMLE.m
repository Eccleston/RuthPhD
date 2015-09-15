function [ml,aic,bic,dic] = ProcessMLE(fname)
%% Analyse the maximum likelihood scores from multichain Filzbach runs

%fname = 'Z:\MHC\Carter2\BasicCarter_option0_final_out.txt';
fid = fopen(fname);

tline = fgetl(fid);

i = 0;
while ischar(tline)
  if ~isempty(tline)
    if strcmp(tline(1:4),'max_')
      i = i + 1;
      ml = textscan(tline,'max_likelihood\t%f'); ml = ml{1};
      tline = fgetl(fid);
      aic = textscan(tline,'AIC\t%f'); aic = aic{1};
      tline = fgetl(fid);
      bic = textscan(tline,'BIC\t%f'); bic = bic{1};
      tline = fgetl(fid);
      dic = textscan(tline,'DIC\t%f'); dic = dic{1};
      continue;
    end
  end
  tline = fgetl(fid);
end

fclose(fid);

return