SSL_none = struct('datafile','SSL_surf_none_data',...
  'sweep','mysweep1',...
  'title','SSL - No IFN',...
  'locs',[8 8 8 8 8 8 8 8],...
  'outname',[outname '_SSL_none']);

SSL_ifn = struct('datafile','SSL_surf_IFN1_data',...
  'sweep','mysweep2',...
  'title','SSL - 1 ng IFN',...
  'locs',[8 8 7 8 7 8 8 8],...
  'outname',[outname '_SSL_ifn']);

ASN_none = struct('datafile','ASN_surf_none_data',...
  'sweep','mysweep3',...
  'title','ASN - No IFN',...
  'locs',[8 8 8 8 8 8 8 8],...
  'outname',[outname '_ASN_none']);

ASN_ifn = struct('datafile','ASN_surf_IFN1_data1',...
  'sweep','mysweep4',...
  'title','ASN - 1 ng IFN',...
  'locs',[8 8 7 8 8 7 7 7],...
  'outname',[outname '_ASN_ifn']);
