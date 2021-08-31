% thesepaths=['/broad/hptmp/lbinan/jeffmicroglia/Analysed/Run1/slice1side2/slice1side2/'];
upperCPN=[18,19,20,21,22,23];
L4stellate=[24,25,26];
deepcpn=[27,28,29];
scpn=[30,31,32,33,34,35,36,37];
chrpn=[38,39,40,41];
l6b=[43,44];
neurons=[45,46];
astrocytes=[47];
endothelium=[48,49,50];
oligo=[51:62];
interneuron=[63,64];
BAM=[65:67];
for alpha=1:size(thesepaths,1)
    mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
    mkdir(savepath);
    tablepath=fullfile(mypath, 'analysis/');
    if alpha==1
        barcodes=table2array(readtable(fullfile(tablepath,strcat(mypath(end-11:end-1),'allcellcounts.csv'))));
        barcodes=[ones(size(barcodes,1),1),barcodes];
    else
        newbarcode=table2array(readtable(fullfile(tablepath,strcat(mypath(end-11:end-1),'allcellcounts.csv'))));
        newbarcode=[alpha*ones(size(newbarcode,1),1),newbarcode];
        barcodes=[barcodes;newbarcode];
    end
end%compile barcodes

% mypath='\\helium\broad_hptmp\lbinan\jeffmicroglia\Analysed\Run1\slice2side2\slice2side2\';
DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
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
% 
modscore=zeros(size(zscored,1),12);
zscoredgenes=zscored(:,5:end);
% %% do it on log(raw)
% zscoredgenes=log(1+o(:,5:end));
% myones=ones(size(o,1),1);
% for i=1:size(modscore,1)
%     for j=1:size(upperCPN,2)
%     modscore(i,1)=modscore(i,1)+zscoredgenes(i,upperCPN(1,j));
%     end
%     for j=1:size(L4stellate,2)
%     modscore(i,2)=modscore(i,2)+zscoredgenes(i,L4stellate(1,j));
%     end
%     for j=1:size(deepcpn,2)
%     modscore(i,3)=modscore(i,3)+zscoredgenes(i,deepcpn(1,j));
%     end
%     for j=1:size(scpn,2)
%     modscore(i,4)=modscore(i,4)+zscoredgenes(i,scpn(1,j));
%     end
%     for j=1:size(chrpn,2)
%     modscore(i,5)=modscore(i,5)+zscoredgenes(i,chrpn(1,j));
%     end
%     for j=1:size(l6b,2)
%     modscore(i,6)=modscore(i,6)+zscoredgenes(i,l6b(1,j));
%     end
%     for j=1:size(neurons,2)
%     modscore(i,7)=modscore(i,7)+zscoredgenes(i,neurons(1,j));
%     end
%     for j=1:size(astrocytes,2)
%     modscore(i,8)=modscore(i,8)+zscoredgenes(i,astrocytes(1,j));
%     end
%     for j=1:size(endothelium,2)
%     modscore(i,9)=modscore(i,9)+zscoredgenes(i,endothelium(1,j));
%     end    
%     for j=1:size(oligo,2)
%     modscore(i,10)=modscore(i,10)+zscoredgenes(i,oligo(1,j));
%     end    
%     for j=1:size(interneuron,2)
%     modscore(i,11)=modscore(i,11)+zscoredgenes(i,interneuron(1,j));
%     end    
%     for j=1:size(BAM,2)
%     modscore(i,12)=modscore(i,12)+zscoredgenes(i,BAM(1,j));
%     end    
% end
% totalscore=sum(modscore);
% for i=1:size(modscore,2)
% modscore1(:,i)=modscore(:,i)/(1+totalscore(1,i));
% end
% 
% % writematrix(modscore,fullfile(thesepaths(1,1:end-25),strcat('scoredallcells.csv')));
% threshupperCPN=prctile(modscore(:,1),90);
% threshL4stellate=prctile(modscore(:,2),95);
% threshdeepcpn=prctile(modscore(:,3),96);
% threshscpn=prctile(modscore(:,4),96);
% threshchrpn=prctile(modscore(:,5),94);
% threshl6b=prctile(modscore(:,6),99.5);
% threshneurons=prctile(modscore(:,7),70);
% threshastrocytes=prctile(modscore(:,8),70);
% threshendothelium=prctile(modscore(:,9),90);
% thresholigo=prctile(modscore(:,10),70);
% threshinterneuron=prctile(modscore(:,11),90);
% threshBAM=prctile(modscore(:,12),80);
% clear scored
% scored(:,1:4)=o(:,1:4);
% 
% scored(:,5:size(mytable,2))=zscoredgenes(:,:);
% myones=ones(size(scored,1),12+size(mytable,2));
% findmax=zeros(size(scored,1),12);
% maxcelltype=zeros(size(scored,1),12);
% maxneurontype=zeros(size(scored,1),12);
% withoutcelltypes=modscore;
% justcelltypes=modscore;
% withoutcelltypes(:,[7:12])=zeros(size(withoutcelltypes,1),6);
% justcelltypes(:,1:6)=zeros(size(justcelltypes,1),6);
% for i=1:size(scored,1)
% findmax(i,:)=modscore(i,:)==max(modscore(i,:));
% maxcelltype(i,:)=justcelltypes(i,:)==max(justcelltypes(i,[7:12]));
% maxneurontype(i,:)=withoutcelltypes(i,:)==max(withoutcelltypes(i,[1:6]));
% end
% findmax=maxcelltype+maxneurontype;
% myones(:,end-11:end)=findmax;
% myones=logical(myones);
% best=modscore;
% for i=1:size(myones,1)
%     best(i,:)=logical(findmax(i,end-11:end)).*modscore(i,:);
%     if sum(best(i,1:6))<2
%         best(i,end-11:end-6)=zeros(1,6);
%     end
%     if sum(best(i,7:12))<2
%         best(i,end-5:end)=zeros(1,6);
%     end
% end
% scored(:,size(mytable,2)+1:size(mytable,2)+12)=best;
% % scored(:,1+size(mytable,2))=myones.*(modscore(:,1)>threshupperCPN);
% % scored(:,2+size(mytable,2))=myones.*(modscore(:,2)>threshL4stellate);
% % scored(:,3+size(mytable,2))=myones.*(modscore(:,3)>threshdeepcpn);
% % scored(:,4+size(mytable,2))=myones.*(modscore(:,4)>threshscpn);
% % scored(:,5+size(mytable,2))=myones.*(modscore(:,5)>threshchrpn);
% % scored(:,6+size(mytable,2))=myones.*(modscore(:,6)>threshl6b);
% % scored(:,7+size(mytable,2))=myones.*(modscore(:,7)>threshneurons);
% % scored(:,8+size(mytable,2))=myones.*(modscore(:,8)>threshastrocytes);
% % scored(:,9+size(mytable,2))=myones.*(modscore(:,9)>threshendothelium);
% % scored(:,10+size(mytable,2))=myones.*(modscore(:,10)>thresholigo);
% % scored(:,11+size(mytable,2))=myones.*(modscore(:,11)>threshinterneuron);
% % scored(:,12+size(mytable,2))=myones.*(modscore(:,12)>threshBAM);
% r=scored;
% r(:,7+size(mytable,2))=0;
% total=sum(r(:,1+size(mytable,2):12+size(mytable,2)),2);
% total=sum(scored(:,1+size(mytable,2):12+size(mytable,2)),2);
% nonattributed=size(total(total(:,1)==0,1));
% attributedseveraltimes=size(total(total(:,1)>1,1));
% for alpha=1:size(thesepaths,1)
%     mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
%     DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
%     zeroImage=zeros(size(DAPI));dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
% 
%     for k=1:12
%         myImage=zeroImage;
% %         for i=1:size(modscore,1)
%             nonzerovalues=fix(zscored(:,4)).*(scored(:,end-12+k)>1);
%             thisregion=zscored(:,1)==alpha;
%             nonzerovalues=nonzerovalues(thisregion);
%             myImage(nonzerovalues(nonzerovalues>0))=1;
% %         end
%         if k==1
%             thisname='uppercpn';
%         end
%         if k==2
%             thisname='l4stellate';
%         end
%         if k==3
%             thisname='deepcpn';
%         end
%         if k==4
%             thisname='scpn';
%         end
%         if k==5
%             thisname='chtpn';
%         end
%         if k==6
%         thisname='l6b';
%         end
%         if k==7
%             thisname='neurons';
%         end
%         if k==8
%             thisname='astrocytes';
%         end
%          if k==9
%             thisname='endothelium';
%         end
%          if k==10
%             thisname='oligo';
%         end        
%          if k==11
%             thisname='interneuron';
%          end   
%         if k==12
%             thisname='BAM';
%         end 
%         maphomeostatic2=cat(3,reduceImage(dapi),255*reduceImage(imdilate(myImage,strel('disk', 20))),reduceImage(zeroImage));
%     imwrite(maphomeostatic2,fullfile(savepath,strcat('Run',num2str(beta),mypath(end-11:end-1),'maplog_',thisname,'.png')));
%     %figure, imshow(maphomeostatic2);
%     end
% end

%% do it again on zscored
zscoredgenes=zscored(:,5:end);
myones=ones(size(o,1),1);
for i=1:size(modscore,1)
    for j=1:size(upperCPN,2)
    modscore(i,1)=modscore(i,1)+zscoredgenes(i,upperCPN(1,j));
    end
    for j=1:size(L4stellate,2)
    modscore(i,2)=modscore(i,2)+zscoredgenes(i,L4stellate(1,j));
    end
    for j=1:size(deepcpn,2)
    modscore(i,3)=modscore(i,3)+zscoredgenes(i,deepcpn(1,j));
    end
    for j=1:size(scpn,2)
    modscore(i,4)=modscore(i,4)+zscoredgenes(i,scpn(1,j));
    end
    for j=1:size(chrpn,2)
    modscore(i,5)=modscore(i,5)+zscoredgenes(i,chrpn(1,j));
    end
    for j=1:size(l6b,2)
    modscore(i,6)=modscore(i,6)+zscoredgenes(i,l6b(1,j));
    end
    for j=1:size(neurons,2)
    modscore(i,7)=modscore(i,7)+zscoredgenes(i,neurons(1,j));
    end
    for j=1:size(astrocytes,2)
    modscore(i,8)=modscore(i,8)+zscoredgenes(i,astrocytes(1,j));
    end
    for j=1:size(endothelium,2)
    modscore(i,9)=modscore(i,9)+zscoredgenes(i,endothelium(1,j));
    end    
    for j=1:size(oligo,2)
    modscore(i,10)=modscore(i,10)+zscoredgenes(i,oligo(1,j));
    end    
    for j=1:size(interneuron,2)
    modscore(i,11)=modscore(i,11)+zscoredgenes(i,interneuron(1,j));
    end    
    for j=1:size(BAM,2)
    modscore(i,12)=modscore(i,12)+zscoredgenes(i,BAM(1,j));
    end    
end
totalscore=sum(modscore);
for i=1:size(modscore,2)
modscore(:,i)=modscore(:,i)/(1+totalscore(1,i));
end
saveallscores=modscore;
% writematrix(modscore,fullfile(thesepaths(1,1:end-25),strcat('scoredallcells.csv')));
clear scored
scored(:,1:4)=o(:,1:4);
% writematrix(modscore,fullfile(thesepaths(1,1:end-25),strcat('scoredallcells.csv')));
threshupperCPN=prctile(modscore(:,1),90);
threshL4stellate=prctile(modscore(:,2),95);
threshdeepcpn=prctile(modscore(:,3),96);
threshscpn=prctile(modscore(:,4),96);
threshchrpn=prctile(modscore(:,5),94);
threshl6b=prctile(modscore(:,6),99.5);
threshneurons=prctile(modscore(:,7),70);
threshastrocytes=prctile(modscore(:,8),70);
threshendothelium=prctile(modscore(:,9),90);
thresholigo=prctile(modscore(:,10),70);
threshinterneuron=prctile(modscore(:,11),90);
threshBAM=prctile(modscore(:,12),80);



scored(:,5:size(mytable,2))=zscoredgenes(:,:);
myones=ones(size(scored,1),12+size(mytable,2));
findmax=zeros(size(scored,1),12);
maxcelltype=zeros(size(scored,1),12);
maxneurontype=zeros(size(scored,1),12);
withoutcelltypes=modscore;
justcelltypes=modscore;
withoutcelltypes(:,[7:12])=zeros(size(withoutcelltypes,1),6);
justcelltypes(:,1:6)=zeros(size(justcelltypes,1),6);
for i=1:size(scored,1)
findmax(i,:)=modscore(i,:)==max(modscore(i,:));
maxcelltype(i,:)=justcelltypes(i,:)==max(justcelltypes(i,[7:12]));
maxneurontype(i,:)=withoutcelltypes(i,:)==max(withoutcelltypes(i,[1:6]));
end
findmax=maxcelltype+maxneurontype;
myones(:,end-11:end)=findmax;
myones=logical(myones);
best=modscore;
for i=1:size(myones,1)
    best(i,:)=logical(findmax(i,end-11:end)).*modscore(i,:);
    if sum(best(i,1:6))<2
        best(i,end-11:end-6)=zeros(1,6);
    end
    if sum(best(i,7:12))<2
        best(i,end-5:end)=zeros(1,6);
    end
end
scored(:,size(mytable,2)+1:size(mytable,2)+12)=best;
% scored(:,1+size(mytable,2))=myones.*(modscore(:,1)>threshupperCPN);
% scored(:,2+size(mytable,2))=myones.*(modscore(:,2)>threshL4stellate);
% scored(:,3+size(mytable,2))=myones.*(modscore(:,3)>threshdeepcpn);
% scored(:,4+size(mytable,2))=myones.*(modscore(:,4)>threshscpn);
% scored(:,5+size(mytable,2))=myones.*(modscore(:,5)>threshchrpn);
% scored(:,6+size(mytable,2))=myones.*(modscore(:,6)>threshl6b);
% scored(:,7+size(mytable,2))=myones.*(modscore(:,7)>threshneurons);
% scored(:,8+size(mytable,2))=myones.*(modscore(:,8)>threshastrocytes);
% scored(:,9+size(mytable,2))=myones.*(modscore(:,9)>threshendothelium);
% scored(:,10+size(mytable,2))=myones.*(modscore(:,10)>thresholigo);
% scored(:,11+size(mytable,2))=myones.*(modscore(:,11)>threshinterneuron);
% scored(:,12+size(mytable,2))=myones.*(modscore(:,12)>threshBAM);
r=scored;
r(:,7+size(mytable,2))=0;
total=sum(r(:,1+size(mytable,2):12+size(mytable,2)),2);
total=sum(scored(:,1+size(mytable,2):12+size(mytable,2)),2);
nonattributed=size(total(total(:,1)==0,1));
attributedseveraltimes=size(total(total(:,1)>1,1));

saveallscores=[saveallscores,best,zscored(:,4),zscored(:,1)];

writematrix(saveallscores,fullfile(thesepaths(1,1:end-25),strcat(mypath(end-11:end-1),'NEWWallcellsscores.csv')));

for alpha=1:size(thesepaths,1)
% for alpha=2:2 
    mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
    DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
    zeroImage=zeros(size(DAPI));
    dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
    mydapi=reduceImage(imbinarize(DAPI,'adaptive','sensitivity',0.3));
    for k=1:12
        myImage=zeroImage;
%         for i=1:size(modscore,1)
            nonzerovalues=fix(zscored(:,4)).*(scored(:,end-12+k)>1);
            thisregion=zscored(:,1)==alpha;
            nonzerovalues=nonzerovalues(thisregion);
            myImage(nonzerovalues(nonzerovalues>0))=1;
%         end
        if k==1
            thisname='uppercpn';
            IMuppercpn=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMuppercpn,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));
        end
        if k==2
            thisname='l4stellate';
            IMl4stellate=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMl4stellate,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));
        end
        if k==3
            thisname='deepcpn';
            IMdeepcpn=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMdeepcpn,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end
        if k==4
            thisname='scpn';
            IMscpn=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMscpn,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));
        end
        if k==5
            thisname='chtpn';
            IMchtpn=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMchtpn,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end
        if k==6
            thisname='l6b';
            IMl6b=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMl6b,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end
        if k==7
            thisname='neurons';
        end
        if k==8
            thisname='astrocytes';
            IMastrocytes=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMastrocytes,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end
         if k==9
            thisname='endothelium';
            IMendothelium=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMendothelium,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end
         if k==10
            thisname='oligo';
            IMoligos=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMoligos,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));

        end        
         if k==11
            thisname='interneuron';
            IMinterneurons=reduceImage(imdilate(myImage,strel('disk', 30)));
            imwrite(IMinterneurons,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));
         end   
        if k==12
            thisname='BAM';
        end 

    end
     redimage=0*IMuppercpn+0.5451*IMl4stellate+1*IMdeepcpn+1*IMscpn+0.9333*IMchtpn+0.4078*IMl6b+0.128*mydapi;
    greenimage=0.6039*IMuppercpn+0.3961*IMl4stellate+0.1882*IMdeepcpn+0.7255*IMscpn+0.1333*IMchtpn+0.651*IMl6b+0.128*mydapi;
    blueimage=0.8039*IMuppercpn+0.0314*IMl4stellate+0.1882*IMdeepcpn+0.0588*IMscpn+0.5451*IMchtpn+0.1176*IMl6b+0.128*mydapi;
    thisImage=cat(3,redimage, greenimage, blueimage);
    imwrite(thisImage,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/neuronmaps/',strcat('Run',num2str(beta),mypath(end-11:end-1),'mapzscored_',thisname,'.png')));
   
end

