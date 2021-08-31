%% do it on zscore

for alpha=1:size(thesepaths,1)
% for alpha=2:2

    mypath=thesepaths(alpha,:)
    savepath=fullfile(mypath, 'analysis/newmosaics/');
    tablepath=fullfile(mypath, 'analysis/');
    mkdir(savepath)
    barcodes=table2array(readtable(fullfile(tablepath,strcat(mypath(end-11:end-1),'microgliacellcountsV2.csv'))));
    disp(num2str(size(barcodes,1)));
    realMg=tdfread(fullfile(mainfolderpath,'kwanhocalls',strcat(thissample(alpha,:),'_microglia_major+PNsubtype.tsv')));
    disp(num2str(sum(realMg.Mg)));
    barcodes=barcodes(logical(realMg.Mg),:);
    barcodes=[alpha*ones(size(barcodes,1),1),barcodes];
    allcells=table2array(readtable(fullfile(strcat('/broad/hptmp/lbinan/jeffmicroglia/Analysed/Run',num2str(beta)),'allcellsscores.csv')));
    allcells=allcells(allcells(:,26)==alpha,:);
    %compile barcodes
DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
% this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
% zeroImage=zeros(this.Height, this.Width);
mytable=barcodes;
this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
zeroImage=zeros(this.Height, this.Width);
for i=1:size(mytable,1)%size normalize
    mytable(i,5:end)=10000*mytable(i,5:end)/(mytable(i,2)*mytable(i,3));
end

o=mytable(mytable(:,2)~=0,:);
for i=1:size(o,1)%totla transcript normalize
    o(i,5:end)=o(i,5:end)/sum(o(i,5:end))*100;
end
zscored=o;
for i=5:size(zscored,2)%zscore
    mymean=mean(zscored(:,i));
    mystd=std(zscored(:,i));
    zscored(:,i)=(zscored(:,i)-mymean)/mystd;
end
  poolzscores=[poolzscores;zscored];
    
    
    
    
zscoredalpha=zscored(zscored(:,1)==alpha,:);
modscore=zeros(size(zscoredalpha,1),8);
oalpha=o(o(:,1)==alpha,:);
zscoredalphagenes=zscoredalpha(:,5:end);
mytablealpha=mytable(mytable(:,1)==alpha,:);
myones=ones(size(oalpha,1),1);
for i=1:size(modscore,1)
    modscore(i,1)=2+zscoredalphagenes(i,1)+zscoredalphagenes(i,2);
    modscore(i,2)=(zscoredalphagenes(i,3)+zscoredalphagenes(i,97)+zscoredalphagenes(i,98)+zscoredalphagenes(i,102)+zscoredalphagenes(i,106)+zscoredalphagenes(i,107)+zscoredalphagenes(i,108))/modscore(i,1)/7;%homeostatic1
    modscore(i,3)=(zscoredalphagenes(i,4)+zscoredalphagenes(i,5)+zscoredalphagenes(i,6)+zscoredalphagenes(i,7)+zscoredalphagenes(i,8)+zscoredalphagenes(i,9)+zscoredalphagenes(i,109)+zscoredalphagenes(i,110)+zscoredalphagenes(i,111)+zscoredalphagenes(i,112))/modscore(i,1)/10;%homeostatic2
    modscore(i,4)=(zscoredalphagenes(i,13)+zscoredalphagenes(i,14))/modscore(i,1)/2;%innate immune
    modscore(i,5)=(zscoredalphagenes(i,11)+zscoredalphagenes(i,12)+zscoredalphagenes(i,100)+zscoredalphagenes(i,101)+zscoredalphagenes(i,105))/modscore(i,1)/5;%inflammatory
    modscore(i,6)=(zscoredalphagenes(i,99)+zscoredalphagenes(i,10))/modscore(i,1)/2;%ApoeHigh
    modscore(i,7)=(zscoredalphagenes(i,103)+zscoredalphagenes(i,114)+zscoredalphagenes(i,115))/modscore(i,1)/3;%Ccr1 high
    modscore(i,8)=(zscoredalphagenes(i,15)+zscoredalphagenes(i,16)+zscoredalphagenes(i,25))/modscore(i,1)/3;%proliferative
end
modscore(:,9)=zscoredalpha(:,1);
% modscore(zscoredalphagenes(:,97)<0,2)=0;modscore(zscoredalphagenes(:,102)<0,2)=0;
% modscore(zscoredalphagenes(:,6)<0,3)=0;modscore(zscoredalphagenes(:,109)<0,3)=0;modscore(zscoredalphagenes(:,111)<0,3)=0;
saveallscores=modscore;
modscore(zscoredalphagenes(:,99)<1,6)=0;
modscore(zscoredalphagenes(:,103)<0.5,7)=0;
modscore(zscoredalphagenes(:,114)<-0.1,7)=0;
modscore(zscoredalphagenes(:,115)<-0.1,7)=0;
modscore(zscoredalphagenes(:,12)<-0.2,5)=0;
modscore(zscoredalphagenes(i,10)<-0.1,6)=0;
modscore(zscoredalphagenes(:,100)<0.1,5)=0;
modscore(zscoredalphagenes(:,101)<0.1,5)=0;

modscore(zscoredalphagenes(:,13)<-0.2,4)=0;
modscore(zscoredalphagenes(:,14)<-0.2,4)=0;
modscore(zscoredalphagenes(:,15)<0,8)=0;
modscore(zscoredalphagenes(:,16)<0,8)=0;
modscore(zscoredalphagenes(:,17)<0,8)=0;

% microglia_thresh=prctile(modscore(:,1),10);
% oalpha=oalpha(modscore(:,1)>microglia_thresh,:);
% saveallscores=saveallscores(modscore(:,1)>microglia_thresh,:);
% zscoredalpha=zscoredalpha(modscore(:,1)>microglia_thresh,:);
% zscoredalphagenes=zscoredalphagenes(modscore(:,1)>microglia_thresh,:);
% mytablealpha=mytablealpha(modscore(:,1)>microglia_thresh,:);
% modscore=modscore(modscore(:,1)>microglia_thresh,:);

for a=1:size(modscore,1)
    confidence(a,:)=modscore(a,:)/sum(modscore(a,:)+min(modscore,[],'all'));
end
writematrix(confidence,fullfile(thesepaths(1,1:end-25),strcat(mypath(end-11:end-1),'confidencemynewscores.csv')));

for i=2:9
    modscore(modscore(:,i)<prctile(modscore(:,i),70),i)=-1;
end
% writematrix(modscore,fullfile(thesepaths(1,1:end-25),strcat('scoredallcells.csv')));

clear scored
scored(:,1:4)=oalpha(:,1:4);

scored(:,5:size(mytablealpha,2))=zscoredalphagenes(:,:);
myones=ones(size(scored,1),7+size(mytablealpha,2));
findmax=zeros(size(scored,1),9);

for i=1:size(scored,1)
findmax(i,:)=modscore(i,:)==max(modscore(i,2:8));
end
findmax=findmax(:,2:end-1);
myones(:,end-6:end)=findmax;
myones=logical(myones);
best=modscore;
for i=1:size(myones,1)
    best(i,end-7:end-1)=logical(findmax(i,end-6:end)).*modscore(i,2:8);
    if sum(best(i,1:end))<0.005
        best(i,end-6:end)=zeros(1,7);
    end
end
for i=1:size(best,1)
    if sum(best(i,2:8))==0
        if saveallscores(i,2)>saveallscores(i,3)
            best(i,2)=saveallscores(i,2);
        else
            best(i,3)=saveallscores(i,3);
        end
    end
end
saveallscores=[saveallscores,best,zscoredalpha(:,4)];

scored(:,size(mytablealpha,2)+1:size(mytablealpha,2)+7)=best(:,2:8);
writematrix(saveallscores,fullfile(thesepaths(1,1:end-25),strcat(mypath(end-11:end-1),'myforcedscores.csv')));
disp('starting plots')

    mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
    DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
    zeroImage=zeros(size(DAPI));dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
    homeo1neighbors=[];
    homeo2neighbors=[];
    immuneneighbors=[];
    inflammatoryneighbor=[];
    apeoneighbors=[];
    Ccr1neighbors=[];
    prolneighbors=[];
    
    for k=1:7
        myImage=zeroImage;
        summed=reduceImage(zeroImage);
%         for i=1:size(modscore,1)
            nonzerovalues=fix(zscoredalpha(:,4)).*(scored(:,end-7+k)>0);
            thisregion=zscoredalpha(:,1)==alpha;
            nonzerovalues=nonzerovalues(thisregion);
            myImage(nonzerovalues(nonzerovalues>0))=1;
            nonzerovalues=nonzerovalues(find(nonzerovalues));
            neighbors=zeros(size(nonzerovalues,1),12);
            
            for thisMG=1:size(nonzerovalues,1)
                [row,column]=ind2sub(size(myImage),nonzerovalues(thisMG));
                distances=zeros(size(allcells,1),1);
                for thiscell=1:size(allcells,1)
                    [cellrow,cellcolumn]=ind2sub(size(myImage),allcells(thiscell,25));
                    distances(thiscell,1)=sqrt((cellrow-row)*(cellrow-row)+(cellcolumn-column)*(cellcolumn-column));
                end
                distances(distances>230)=0;
                if sum(distances)>0
                celltokeep=zeros(size(distances,1),1);
                celltokeep(find(distances))=1;
                cellIDs=allcells(:,13:24);
                mycellIDs=zeros(size(cellIDs));
                mycellIDs(find(cellIDs))=1;
                neighbors(thisMG,:)=sum(mycellIDs(logical(celltokeep),:),1);
                end
            end
    if k==1
            thisname='homeostatic1';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==2
            thisname='homeostatic2';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==3
            thisname='innate_immune';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==4
            thisname='inflammatory';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==5
            thisname='ApoeHigh';
            ApoeHigh=reduceImage(imdilate(myImage,strel('disk', 30)));
           writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==6
        thisname='Ccr1High';
       writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));
        end
        if k==7
            thisname='proliferative';
        writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.csv')));

        end
%         if k==1
%             thisname='homeostatic1';
%             homeo1=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+homeo1;
%             imwrite(homeo1,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
%         if k==2
%             thisname='homeostatic2';
%             homeo2=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+homeo2;
%             imwrite(homeo2,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
%         if k==3
%             thisname='innate_immune';
%             innate_immune=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+innate_immune;
%             imwrite(innate_immune,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
%         if k==4
%             thisname='inflammatory';
%             inflammatory=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+inflammatory;
%             imwrite(inflammatory,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
%         if k==5
%             thisname='ApoeHigh';
%             ApoeHigh=reduceImage(imdilate(myImage,strel('disk', 30)));
%             imwrite(ApoeHigh,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
%         if k==6
%         thisname='Ccr1High';
%         Ccr1High=reduceImage(imdilate(myImage,strel('disk', 30)));
%         summed=summed+Ccr1High;
%         imwrite(Ccr1High,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         
%         end
%         if k==7
%             thisname='proliferative';
%         proliferative=reduceImage(imdilate(myImage,strel('disk', 30)));
%          summed=summed+proliferative;
%          imwrite(proliferative,fullfile('/broad/hptmp/lbinan/jeffmicroglia/prettymapsFORCEDwihKMg/',strcat('run', num2str(beta),mypath(end-11:end-1),thisname,'.png')));
% 
%         end
    end
%     onlyOnes=ones(size(summed));
%         onlyOnes(summed>1)=0;
%         onlyOnes=imbinarize(onlyOnes);
%         redimage=homeo1+0.8*homeo2+0.1216*ApoeHigh+0.6510*innate_immune+0.4588*Ccr1High+0.4*inflammatory+0.9686*proliferative;
%         greenimage=0.8*homeo1+0*homeo2+0.4706*ApoeHigh+0.8078*innate_immune+0.4392*Ccr1High+0.651*inflammatory+0.5059*proliferative;
%         blueimage=0.4*homeo1+0*homeo2+0.7059*ApoeHigh+0.8902*innate_immune+0.7020*Ccr1High+0.1176*inflammatory+0.7490*proliferative;
%         redimage=redimage.*onlyOnes;
%         greenimage=greenimage.*onlyOnes;
%         blueimage=blueimage.*onlyOnes;
%         thisImage=cat(3,redimage, greenimage, blueimage);
%         imwrite(thisImage,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/heximages/threshold/',strcat('run', num2str(beta),mypath(end-11:end-1),'new_70.png')));


end
% save('poolzscores.mat','poolzscores');