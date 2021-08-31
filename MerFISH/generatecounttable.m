function tizdone=generatecounttable(mypath)
%mypath='\\helium\broad_hptmp\lbinan\jeffmicroglia\Analysed\Run1\slice2side2\slice2side2\';
savepath=fullfile(mypath, 'analysis');
mkdir(savepath);
if isfile(fullfile(savepath,'cellmask.mat'))
    load(fullfile(savepath,'cellmask.mat'))
else
    if strcmp(mypath(1,44:47),'Run2')
        cellmask=findcells(mypath);
    elseif strcmp(mypath(1,44:47),'Run1')
        cellmask=findcellsforRun1(mypath);
    end
    save(fullfile(savepath,'cellmask.mat'),'cellmask');
end
stats=regionprops(cellmask,'area','PixelIdxList','centroid');
zeroImage=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
outputpath=fullfile(mypath,'Merfish_mosaics');
Zisbarcodes=table2array(readtable(fullfile(savepath,'allbarcodeswithcellID.csv')));
zeroImage=zeros(size(zeroImage));
countable=zeros(size(stats,1),115+3);
initline=countable(1,:);
for i=1:size(stats,1)
    for z=1:7
        mytable=Zisbarcodes(Zisbarcodes(:,2)==z,:);
        mytable=mytable(mytable(:,3)==i,:);
        gene1=size(mytable(mytable(:,1)==1,:),1);
        gene2=size(mytable(mytable(:,1)==2,:),1);
        microglialgenes=gene1+gene2;
        if microglialgenes>=6
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

writematrix(countable,fullfile(savepath,strcat(mypath(end-11:end-1),'cellcounts.csv')));
%% all cells

if isfile(fullfile(savepath,'allcellmask.mat'))
    load(fullfile(savepath,'allcellmask.mat'))
else
    if strcmp(mypath(1,44:47),'Run2')
        cellmask=findallcells(mypath,genelist);
    elseif strcmp(mypath(1,44:47),'Run1')
        cellmask=findallcellsforRun1(mypath,genelist));
    end
    save(fullfile(savepath,'allcellmask.mat'),'cellmask');
end
Zisbarcodes=table2array(readtable(fullfile(savepath,'allcellsbarcodes.csv')));
zeroImage=zeros(size(zeroImage));
countable=zeros(size(stats,1),115+3);
initline=countable(1,:);

for i=1:size(stats,1)
    for z=1:7
        mytable=Zisbarcodes(Zisbarcodes(:,2)==z,:);
        mytable=mytable(mytable(:,3)==i,:);
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

writematrix(countable,fullfile(savepath,strcat(mypath(end-11:end-1),'allcellcounts.csv')));
tizdone='done';