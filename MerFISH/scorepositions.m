function mytable=scorepositions(thissample,zeroImage,cellimage,mypath)
% mypath='\\helium\broad_hptmp\lbinan\jeffmicroglia\Analysed\Run2\slice3side2\slice3side2\';
load(fullfile(mypath,'positions2keep.mat'));
% load('\\helium\broad_hptmp\lbinan\jeffmicroglia\Analysed\angletable.mat');
load('/broad/hptmp/lbinan/jeffmicroglia/Analysed/angletable.mat');
myangle=angletable(thissample,1);
myimage=zeroImage;
myimage(keeppositions)=1;
myimage=imrotate(myimage,myangle);
cellimage=imrotate(cellimage,myangle);
mymin=0;
% for j=i:size(myimage,2)
%     if myimage(fix(size(myimage,1)/2),j)~=0
%         mymin=j;
%         break
%     end
% end
% mymax=0;
mymax=0;
mynewmax=0;
for i=1:size(myimage,1)
for j=1:size(myimage,2)-1
    if myimage(i,size(myimage,2)-j)~=0
        mynewmax=size(myimage,2)-j;
        break
    end
end
mymax=max(mymax,mynewmax);
end
mypositions=find(cellimage);
mytable=[];
for c=1:size(mypositions,1)
    cellposition=mypositions(c,1);
    [row,col] = ind2sub(size(cellimage),cellposition);
%     myref=0;
%     for y=1:size(cellimage,2)
%         if myimage(row,y)==1
%         myref=y;
%         end
%     end
    thisscore=mymax-col;
    mytable=[mytable;thisscore];
end