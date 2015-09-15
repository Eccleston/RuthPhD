%% Example of generating string

base = fileread('sweep_data.lbs');%('lbs_base.lbs');
datadir = '../Data/';

dat1 = dlmread([datadir 'AllData_150520_ASN_SSL.txt'],',',1,1);
locs = 3:66;
SSL_cyt_none = reshape(dat1(locs,1),8,8);
SSL_surf_none = reshape(dat1(locs,2),8,8);
SSL_cyt_ifn1 = reshape(dat1(locs,3),8,8);
SSL_surf_ifn1 = reshape(dat1(locs,4),8,8);
SSL_cyt_ifn2 = reshape(dat1(locs,5),8,8);
SSL_surf_ifn2 = reshape(dat1(locs,6),8,8);
ASN_cyt_none = reshape(dat1(locs,7),8,8);
ASN_cyt_ifn1 = reshape(dat1(locs,8),8,8);
ASN_cyt_ifn2 = reshape(dat1(locs,9),8,8);

A = [SSL_cyt_none(:,1) ASN_cyt_none(:,1);
     SSL_cyt_none(:,2) ASN_cyt_none(:,2);
     SSL_cyt_none(:,3) ASN_cyt_none(:,3);
     SSL_cyt_none(:,4) ASN_cyt_none(:,4);
     SSL_cyt_none(:,5) ASN_cyt_none(:,5);
     SSL_cyt_none(:,6) ASN_cyt_none(:,6);
     SSL_cyt_none(:,7) ASN_cyt_none(:,7);
     SSL_cyt_none(:,8) ASN_cyt_none(:,8);];
%A = [1.0 1.0;
 % 2.0 1.0;
  %3.0 1.0;
  %4.0 1.0;
  %5.0 1.0;
  %1.0 2.0;
  %2.0 2.0;
  %3.0 2.0;
  %4.0 2.0;
  %5.0 2.0];

str = 'directive sweep mysweep1 = {upreg = [0.0], (in1,in2) = [';
% Start with first entry
str = sprintf('%s(%1.2f,%1.2f)',str,A(1,1),A(1,2));
% Then add a comma and additional entries in a loop
for i = 2:length(A)
  str = sprintf('%s,(%1.2f,%1.2f)',str,A(i,1),A(i,2));
end
% Add the final delimiter
str = strcat(str,'] }');

B = [SSL_cyt_ifn1(:,1) ASN_cyt_ifn1(:,1);
     SSL_cyt_ifn1(:,2) ASN_cyt_ifn1(:,2);
     SSL_cyt_ifn1(:,3) ASN_cyt_ifn1(:,3);
     SSL_cyt_ifn1(:,4) ASN_cyt_ifn1(:,4);
     SSL_cyt_ifn1(:,5) ASN_cyt_ifn1(:,5);
     SSL_cyt_ifn1(:,6) ASN_cyt_ifn1(:,6);
     SSL_cyt_ifn1(:,7) ASN_cyt_ifn1(:,7);
     SSL_cyt_ifn1(:,8) ASN_cyt_ifn1(:,8);];
 

validdata_SSL_surf_ifn1 = ~isnan(SSL_surf_ifn1);
validdata_SSL_cyt_ifn1_1 = ~isnan(SSL_cyt_ifn1);
validdata_ASN_cyt_ifn1_1 = ~isnan(ASN_cyt_ifn1);

validdataAll = validdata_SSL_surf_ifn1 & validdata_SSL_cyt_ifn1_1  & validdata_ASN_cyt_ifn1_1  ;

strB = 'directive sweep mysweep2 = {upreg = [1.0], (in1,in2) = [';
% Start with first entry
strB = sprintf('%s(%1.2f,%1.2f)',strB,B(1,1),B(1,2));
% Then add a comma and additional entries in a loop
for i = 2:length(b)
  strB = sprintf('%s,(%1.2f,%1.2f)',strB,B(i,1),B(i,2));
end
% Add the final delimiter
strB = strcat(strB,'] }');


fid = fopen('temp.lbs','w');
fprintf(fid,'//Auto-generated sweep\r\n%s\r\n\r\n',str);

fprintf(fid,'//Auto-generated sweep\r\n%s\r\n\r\n',strB);

fclose(fid);
system('copy /b temp.lbs+150520_ASN_SSL_sweep.lbs combined_new3.lbs');
%system('copy /b temp.lbs+lbs_base.lbs combined.lbs');

for k=1:8
    for j=1:8
        str4='SSL_';
eval(sprintf('%s%1.2d_%1.2d=%1.2d',str4,SSL_cyt_ifn1(j,k), ASN_cyt_ifn1(j,k),SSL_surf_ifn1(j,k)));
    end
end
str5='Time';
str5=sprintf('%s \t %s%1.2d_%1.2d',str5,'SSL_',B(1,1),B(1,2));
for i = 2:length(A)
  str5 = sprintf('%s \t %s%1.2d_%1.2d',str5,'SSL_',B(i,1),B(i,2));
end
time=24*3600;
T=zeros(1,1); T=time;
D2=horzcat(SSL_surf_ifn1(:,1)',SSL_surf_ifn1(:,2)',SSL_surf_ifn1(:,3)',SSL_surf_ifn1(:,4)',SSL_surf_ifn1(:,5)',SSL_surf_ifn1(:,6)',SSL_surf_ifn1(:,7)',SSL_surf_ifn1(:,8)')
fid4=fopen('test2.txt','w');
fprintf(fid4, str5);
fprintf(fid4, '\n');
fprintf(fid4, '%1.2d \t',D2);


fclose(fid4);

dat2 = dlmread([datadir 'AllData_150520_ASN_SSL_reverse.txt'],',',1,1);
locs = 3:66;
SSL_cyt_none_2 = reshape(dat2(locs,1),8,8);
ASN_surf_none = reshape(dat2(locs,2),8,8);
SSL_cyt_ifn1_2 = reshape(dat2(locs,3),8,8);
ASN_surf_ifn1 = reshape(dat2(locs,4),8,8);
SSL_cyt_ifn2_2 = reshape(dat2(locs,5),8,8);
ASN_surf_ifn2 = reshape(dat2(locs,6),8,8);
ASN_cyt_none_2 = reshape(dat2(locs,7),8,8);
ASN_cyt_ifn1_2 = reshape(dat2(locs,8),8,8);
ASN_cyt_ifn2_2 = reshape(dat2(locs,9),8,8);

C = [SSL_cyt_none_2(:,1) ASN_cyt_none_2(:,1);
     SSL_cyt_none_2(:,2) ASN_cyt_none_2(:,2);
     SSL_cyt_none_2(:,3) ASN_cyt_none_2(:,3);
     SSL_cyt_none_2(:,4) ASN_cyt_none_2(:,4);
     SSL_cyt_none_2(:,5) ASN_cyt_none_2(:,5);
     SSL_cyt_none_2(:,6) ASN_cyt_none_2(:,6);
     SSL_cyt_none_2(:,7) ASN_cyt_none_2(:,7);
     SSL_cyt_none_2(:,8) ASN_cyt_none_2(:,8);];

D = [SSL_cyt_ifn1_2(:,1) ASN_cyt_ifn1_2(:,1);
     SSL_cyt_ifn1_2(:,2) ASN_cyt_ifn1_2(:,2);
     SSL_cyt_ifn1_2(:,3) ASN_cyt_ifn1_2(:,3);
     SSL_cyt_ifn1_2(:,4) ASN_cyt_ifn1_2(:,4);
     SSL_cyt_ifn1_2(:,5) ASN_cyt_ifn1_2(:,5);
     SSL_cyt_ifn1_2(:,6) ASN_cyt_ifn1_2(:,6);
     SSL_cyt_ifn1_2(:,7) ASN_cyt_ifn1_2(:,7);
     SSL_cyt_ifn1_2(:,8) ASN_cyt_ifn1_2(:,8);];


 
base_new = fileread('check11.lbs');
 
str6 = 'directive sweep mysweep3 = {upreg = [0.0], (in1,in2) = [';
% Start with first entry
str6 = sprintf('%s(%1.2f,%1.2f)',str6,C(1,1),C(1,2));
% Then add a comma and additional entries in a loop
for i = 2:length(C)
  str6 = sprintf('%s,(%1.2f,%1.2f)',str6,C(i,1),C(i,2));
end
% Add the final delimiter
str6 = strcat(str6,'] }');


validdata_SSL_surf_ifn1 = ~isnan(SSL_surf_ifn1);
validdata_SSL_cyt_ifn1_1 = ~isnan(SSL_cyt_ifn1);
validdata_ASN_cyt_ifn1_1 = ~isnan(ASN_cyt_ifn1);

validdataAll = validdata_SSL_surf_ifn1 & validdata_SSL_cyt_ifn1_1  & validdata_ASN_cyt_ifn1_1  ;

str7 = 'directive sweep mysweep4 = {upreg = [1.0], (in1,in2) = [';
% Start with first entry
str7 = sprintf('%s(%1.2f,%1.2f)',str7,D(1,1),D(1,2));
% Then add a comma and additional entries in a loop
for i = 2:length(D)
  str7 = sprintf('%s,(%1.2f,%1.2f)',str7,D(i,1),D(i,2));
end
% Add the final delimiter
str7 = strcat(str7,'] }');
fid3 = fopen('temp.lbs','w');
fprintf(fid3,'//Auto-generated sweep\r\n%s\r\n\r\n',str6);

fprintf(fid3,'//Auto-generated sweep\r\n%s\r\n\r\n',str7);

fclose(fid3);
system('copy /b temp.lbs+check11.lbs combined_new_test.lbs');

for k=1:8
    for j=1:8
        str8='ASN_';
eval(sprintf('%s%1.2d_%1.2d=%1.2d',str8,SSL_cyt_ifn1_2(j,k), ASN_cyt_ifn1_2(j,k),ASN_surf_ifn1(j,k)));
    end
end
str9='Time';
str9=sprintf('%s \t %s%1.2d_%1.2d',str9,'ASN_',D(1,1),D(1,2));
for i = 2:length(D)
  str9 = sprintf('%s \t %s%1.2d_%1.2d',str9,'ASN_',D(i,1),D(i,2));
end
time=24*3600;
T=zeros(1,1); T=time;
D2=horzcat(ASN_surf_ifn1(:,1)',ASN_surf_ifn1(:,2)',ASN_surf_ifn1(:,3)',ASN_surf_ifn1(:,4)',ASN_surf_ifn1(:,5)',ASN_surf_ifn1(:,6)',ASN_surf_ifn1(:,7)',ASN_surf_ifn1(:,8)')
fid6=fopen('test4.txt','w');
fprintf(fid6, str9);
fprintf(fid6, '\n');
fprintf(fid6, '%1.2d \t',D2);


fclose(fid6);


for k=1:8
    for j=1:8
        str10='ASN_';
eval(sprintf('%s%1.2d_%1.2d=%1.2d',str10,SSL_cyt_none_2(j,k), ASN_cyt_none_2(j,k),ASN_surf_none(j,k)));
    end
end
str11='Time';
str11=sprintf('%s \t %s%1.2d_%1.2d',str11,'ASN_',C(1,1),C(1,2));
for i = 2:length(C)
  str11 = sprintf('%s \t %s%1.2d_%1.2d',str11,'ASN_',C(i,1),C(i,2));
end
time=24*3600;
T=zeros(1,1); T=time;
D2=horzcat(ASN_surf_none(:,1)',ASN_surf_none(:,2)',ASN_surf_none(:,3)',ASN_surf_none(:,4)',ASN_surf_none(:,5)',ASN_surf_none(:,6)',ASN_surf_none(:,7)',ASN_surf_none(:,8)')
fid7=fopen('test6.txt','w');
fprintf(fid7, str11);
fprintf(fid7, '\n');
fprintf(fid7, '%1.2d \t',D2);


fclose(fid7);



return