function tizdone=generatecounttableallgenes(mypath)
%mypath='\\helium\broad_hptmp\lbinan\jeffmicroglia\Analysed\Run1\slice2side2\slice2side2\';
savepath=fullfile(mypath, 'analysis');
mkdir(savepath);
load(fullfile(savepath,'allcellmask.mat'))
%cellmask=findcells(mypath);
stats=regionprops(cellmask,'area','PixelIdxList','centroid');
this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
zeroImage=zeros(this.Height, this.Width);
outputpath=fullfile(mypath,'Merfish_mosaics');
Zisbarcodes=table2array(readtable(fullfile(savepath,'allcellsbarcodes.csv')));
countable=zeros(size(stats,1),115+3);
initline=countable(1,:);
for i=1:size(stats,1)
    for z=1:7
        mytable=Zisbarcodes(Zisbarcodes(:,2)==z,:);
        mytable=mytable(mytable(:,3)==i,:);
        if size(mytable,1)>0
        gene1=size(mytable(mytable(:,1)==1,:),1);
        gene2=size(mytable(mytable(:,1)==2,:),1);
            countable(i,1)=countable(i,1)+1;%#of z in that cell
            countable(i,2)=stats(i).Area;
            countable(i,3)=mytable(1,4);%position in image (of first found barcode)
            countable(i,4)=countable(i,4)+gene1;
            countable(i,5)=countable(i,5)+gene2;
            for j=3:115
                thiscount=size(mytable(mytable(:,1)==j,:),1);
                countable(i,j+3)=countable(i,j+3)+thiscount;
            end
        end
    end
end

writematrix(countable,fullfile(savepath,strcat(mypath(end-11:end-1),'allcellcounts.csv')));
tizdone='done';