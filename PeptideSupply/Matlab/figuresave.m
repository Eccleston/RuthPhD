function figuresave(fh,fname,dpi)
%% FIGURESAVE(FH,FNAME,DPI)
%    Saves figure in handle FH to several files
if ~(exist('Figures','dir')==7)
  !mkdir Figures
end

set(fh,'PaperPositionMode','auto')
save2pdf(['Figures/' fname '.pdf'],fh,dpi)

%saveas(fh,['Figures/' fname '.fig'],'fig')
%print('-depsc2','-r300',['Figures/' fname])
%print('-dtiff','-r300',['Figures/' fname])
set(fh,'PaperPositionMode','auto')
print('-dpng',['-r' num2str(dpi)],['Figures/' fname '.png'])

return