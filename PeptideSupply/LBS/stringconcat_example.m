%% Example of generating string

base = fileread('lbs_base.lbs');

A = [1.0 1.0;
  2.0 1.0;
  3.0 1.0;
  4.0 1.0;
  5.0 1.0;
  1.0 2.0;
  2.0 2.0;
  3.0 2.0;
  4.0 2.0;
  5.0 2.0];

str = 'directive sweep mysweep = { (in1,in2) = [';
% Start with first entry
str = sprintf('%s(%1.2f,%1.2f)',str,A(1,1),A(1,2));
% Then add a comma and additional entries in a loop
for i = 2:length(A)
  str = sprintf('%s,(%1.2f,%1.2f)',str,A(i,1),A(i,2));
end
% Add the final delimiter
str = strcat(str,'] }');

fid = fopen('temp.lbs','w');
fprintf(fid,'//Auto-generated sweep\r\n%s\r\n\r\n',str);
fclose(fid);

system('copy /b temp.lbs+lbs_base.lbs combined.lbs');

return