function fh = PlotCorrelation(paras,post)
%%
Np = length(paras);
C = corrcoef(post); Cr = reshape(C,Np^2,1);

Nc = 50;
colors = [[ones(Nc+1,1) (0:1/Nc:1)'*[1 1]];
  [(1:-1/Nc:0)'*[1 1] ones(Nc+1,1)]] .^ 0.5;

X = ones(Np,1)*(1:Np); Xr = reshape(X,Np^2,1);
Y = X'; Yr = reshape(Y,Np^2,1);
locs = find(Xr<Yr);

dx = 75;
dy = 60;
cx = 10;
fwidth = dx + Np*25 + 4*cx;
fheight = dy + Np*25;
cwidth = cx/fwidth;
left = dx/fwidth; width = 1-left-4*cwidth;
bottom = dy/fheight; height = 0.95-bottom;

fh = figure;
fh.Position = [700 100 fwidth fheight];
scatter3(Xr(locs),-Yr(locs),Cr(locs),400,Cr(locs),'.')
axis([0 Np+1 -Np-1 0 -1 1 -1 1])
view([0 90])
set(gca,'Position',[left bottom width height])
set(gca,'Xtick',1:Np,'Xticklabel',paras)
set(gca,'Ytick',-Np:-1,'Yticklabel',fliplr(paras))
rotateticklabel(gca,45,0.02);
colormap(colors)
hc = colorbar;
hc.Limits = [-1 1];
hc.Position = [1-4*cwidth bottom cwidth height];

return