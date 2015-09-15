%% Process the inference from LBS
function post = ProcessParameters(f,logs)

[headers,dataIn] = importTSV([f.simDir 'posterior.tsv']);

% Trim fixed parameters
locs = std(dataIn) > 1e-10;
paras0 = headers(locs);
data = dataIn(:,locs);
for i = 1:length(paras0)
  post.(paras0{i}) = data(:,i);
end
paras = paras0(cellfun('isempty',regexp(paras0,'noise*')));
paras = paras(cellfun('isempty',regexp(paras,'iter*')));
Np = length(paras);
[~,~,ib] = intersect(paras,paras0);

%% Plot the parameter posterior distributions

Nc = 5;
Nr = ceil(Np/Nc);
dx = 0.195;
dy = 1/Nr;
left = 0.05;
bottom = 0.16/Nr;
width = 0.14;
height = 0.68/Nr;
f1 = figure;
f1.Position = [600 500 800 175*Nr];
for i = 1:Np
  row = Nr - ceil(i/Nc);
  col = mod(i-1,Nc);
  subplot('position',[left+col*dx bottom+row*dy width height])
  if ismember(paras{i},logs)
    hist(log10(post.(paras{i})),30)
    title(['log_{10}(' paras{i} ')'])
  else
    hist(post.(paras{i}),30)
    title(paras{i})
  end
  box off
end

sd = sum(abs(data(:,3:end)),2);
acc = length(find(diff(sd)>0)) / (size(data,1)-1);

fprintf('MLE: %1.2f \t Accepted rate: %1.2f\n', max(post.likelihood), acc);
if f.saving
  figuresave(f1,[f.outname '_posterior'],300)
end

%% Corrleations
f2 = PlotCorrelation(paras(2:end),data(1:end,ib(2:end)));

if f.saving
  figuresave(f2,[f.outname '_correlation'],300)
end

return