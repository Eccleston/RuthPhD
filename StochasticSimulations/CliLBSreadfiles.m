% Start with a folder and get a list of all subfolders.
% Finds and prints names of all PNG, JPG, and TIF images in 
% that folder and all of its subfolders.
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;

% Define a starting folder.
start_path = fullfile('C:\Users\RuthCharlotte\Documents\MATLAB\All\All\CliLBS\');
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames)
data_all=zeros(9002,10,numberOfFolders);
% Process all tsv files in those folders.
for k = 2 : numberOfFolders-1
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	fprintf('Processing folder %s\n', thisFolder);
	
	% Get tsv files.
	filePattern = sprintf('%s/*.tsv', thisFolder);
	baseFileNames = dir(filePattern);
	numberOftsvFiles = length(baseFileNames);
	% Now we have a list of all files in this folder.
	data = dlmread([thisFolder '\stochastic.tsv'],'\t',1,0);
    data_all(:,:,k-1)=data;
    time = data(:,1);
    SP = data(:,2:10);
    figure(1);plot(time,SP);hold on;
    
    
 	if numberOftsvFiles >= 1
		% Go through all those tsv files.
		for f = 1 : numberOftsvFiles
			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
			fprintf('     Processing tsv file %s\n', fullFileName);
            
		end
	else
		fprintf('     Folder %s has no tsv files in it.\n', thisFolder);
    end
end
rev=zeros(9002,numberOfFolders);
nef=zeros(9002,numberOfFolders);
tat=zeros(9002,numberOfFolders);
gag=zeros(9002,numberOfFolders);
pol=zeros(9002,numberOfFolders);
env=zeros(9002,numberOfFolders);
vif=zeros(9002,numberOfFolders);
vpr=zeros(9002,numberOfFolders);
vpu=zeros(9002,numberOfFolders);

rev(:,:)=data_all(:,2,:);
nef(:,:)=data_all(:,3,:);
tat(:,:)=data_all(:,4,:);
gag(:,:)=data_all(:,5,:);
pol(:,:)=data_all(:,6,:);
env(:,:)=data_all(:,7,:);
vif(:,:)=data_all(:,8,:);
vpr(:,:)=data_all(:,9,:);
vpu(:,:)=data_all(:,10,:);

rm=mean(rev')';
nm=mean(nef')';
tm=mean(tat')';
gm=mean(gag')';
pm=mean(pol')';
em=mean(env')';
vim=mean(vif')';
vpum=mean(vpu')';
vprm=mean(vpr')';
time=time/3600;
figure(2);hold on
plot(time,nm,time,gm,time,pm,time,em,time,vim,time,vprm,time,vpum,'LineWidth',1.5);legend({'Nef','Gag','Pol','Env','Vif','Vpr','Vpu'});xlabel('time (h)');ylabel('cell surface abundance');set(gca,'FontSize',16)

workingdir = 'C:\Users\RuthCharlotte\Documents\MATLAB\All\All\';
std_gag=std(gag')';
var_gag=var(gag')';
std_pol=std(pol')';
var_pol=var(pol')';
std_env=std(env')';
var_env=var(env')';
std_vpr=std(vpr')';
var_vpr=var(vpr')';
std_nef=std(nef')';
var_nef=var(nef')';
std_vif=std(vif')';
var_vif=var(vif')';
std_vpu=std(vpu')';
var_vpu=var(vpu')';

figure(3);hold on
H = shadedErrorBar(time, gm, std_gag, {'-b', 'LineWidth', 2}, 0); title('Average Gag Peptide Presentation, nVirion=5','FontSize',16);xlabel('time (h)','FontSize',16);ylabel('cell surface abundance','FontSize',16)
legend([H.mainLine, H.patch], ...
'\mu', '\sigma', ...
'Location', 'Northwest');
set(gca,'FontSize',16);

figure(4);hold on
H = shadedErrorBar(time, pm, std_pol, {'-b', 'LineWidth', 2}, 0); title('Average Pol Peptide Presentation, nVirion=5','FontSize',16);xlabel('time (h)','FontSize',16);ylabel('cell surface abundance','FontSize',16)
legend([H.mainLine, H.patch], ...
'\mu', '\sigma', ...
'Location', 'Northwest');
set(gca,'FontSize',16)

figure(5);hold on
H = shadedErrorBar(time, vprm, std_vpr, {'-b', 'LineWidth', 2}, 0); title('Average Vpr Peptide Presentation, nVirion=5','FontSize',16);xlabel('time (h)','FontSize',16);ylabel('cell surface abundance','FontSize',16)
legend([H.mainLine, H.patch], ...
'\mu', '\sigma', ...
'Location', 'Northwest');
set(gca,'FontSize',16)
return
