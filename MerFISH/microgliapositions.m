%% do it on zscore
remember=[];
for alpha=1:size(thesepaths,1)
    
    
    mypath=thesepaths(alpha,:)
    savepath=fullfile(mypath, 'analysis/newmosaics/');
    tablepath=fullfile(mypath, 'analysis/');
    mkdir(savepath)
    barcodes=table2array(readtable(fullfile(tablepath,strcat(mypath(end-11:end-1),'microgliacellcountsV2.csv'))));
    barcodes=[alpha*ones(size(barcodes,1),1),barcodes];
    disp(num2str(size(barcodes,1)));
    realMg=tdfread(fullfile(mainfolderpath,'kwanhocalls',strcat(thatsample(alpha,:),'_microglia_major+PNsubtype.tsv')));
    disp(num2str(sum(realMg.Mg)));
    barcodes=barcodes(logical(realMg.Mg),:);
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
    modscore(i,8)=(zscoredalphagenes(i,15)+zscoredalphagenes(i,16)+zscoredalphagenes(i,17))/modscore(i,1)/3;%proliferative
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
%%

disp('starting plots')

    mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
    DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
    zeroImage=zeros(size(DAPI));dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
%                 load(fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/layers/',strcat('Run', num2str(beta),mypath(end-11:end-1),'L1.mat')));
%             load(fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/layers/',strcat('Run', num2str(beta),mypath(end-11:end-1),'L2_4.mat')));
%             load(fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/layers/',strcat('Run', num2str(beta),mypath(end-11:end-1),'L5.mat')));
%             load(fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/layers/',strcat('Run', num2str(beta),mypath(end-11:end-1),'L6.mat')));
%             Mgcounts=zeros(8,4);
%             Mgcounts(8,1)=size(L1,1);
%             Mgcounts(8,2)=size(L2_4,1);
%             Mgcounts(8,3)=size(L5,1);
%             Mgcounts(8,4)=size(L6,1);
            if beta==1
                thissample=alpha;
            elseif beta==2
                %thissample=alpha+2;
                thissample=alpha+4;
            else
                thissample=alpha+7;
            end
            Thomeo1=[];
Thomeo2=[];
TApeoHigh=[];
TCcr1=[];

Timmune=[];
Tinflammatory=[];
Tproliferative=[];
    for k=1:7
        myImage=zeroImage;
        summed=reduceImage(zeroImage);
%         for i=1:size(modscore,1)
            nonzerovalues=fix(zscoredalpha(:,4)).*(scored(:,end-7+k)>0);
            thisregion=zscoredalpha(:,1)==alpha;
            nonzerovalues=nonzerovalues(thisregion);
            myImage(nonzerovalues(nonzerovalues>0))=1;
        if k==1
                thisname='homeostatic1';
                Thomeo1=[Thomeo1;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==2
            thisname='homeostatic2';
            Thomeo2=[Thomeo2;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==3
            thisname='innate_immune';
            Timmune=[Timmune;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==4
            thisname='inflammatory';
            Tinflammatory=[Tinflammatory;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==5
            thisname='ApoeHigh';
            TApeoHigh=[TApeoHigh;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==6
        thisname='Ccr1High';
       TCcr1=[TCcr1;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==7
            thisname='proliferative';
            Tproliferative=[Tproliferative;scorepositions(thissample,zeroImage,myImage,mypath)];
        end

    end
    
    mysize=max([size(Thomeo1,1),size(Thomeo2,1),size(Timmune,1)]);
    mydistances=zeros(mysize,7);
    mydistances(1:size(Thomeo1,1),1)=Thomeo1;
    mydistances(1:size(Thomeo2,1),2)=Thomeo2;
    mydistances(1:size(Timmune,1),3)=Timmune;
    mydistances(1:size(Tinflammatory,1),4)=Tinflammatory;
    mydistances(1:size(TApeoHigh,1),5)=TApeoHigh;
    mydistances(1:size(TCcr1,1),6)=TCcr1;
    mydistances(1:size(Tproliferative,1),7)=Tproliferative;
    remember=[remember;mydistances];
    writematrix(mydistances,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/piaDistKMg/',strcat('Run', num2str(beta),mypath(end-11:end-1),'myMgdistances.csv')));
    end
writematrix(remember,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/',strcat('Run', num2str(beta),'myMgdistances.csv')));
