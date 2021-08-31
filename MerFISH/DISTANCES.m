%% do it on zscore
mainfolderpath=['/broad/hptmp/lbinan/jeffmicroglia/'];
% mainfolderpath=['\\helium\broad_hptmp\lbinan\jeffmicroglia\'];
thissample=['Run2slice2side2'];
mynamethissample=thesenames;
for alpha=1:size(thissample,1)
    mypath=fullfile(mainfolderpath,'Analysed','Run2',thissample(alpha,5:end),thissample(alpha,5:end));
    microgliapositions=table2array(readtable(fullfile(mainfolderpath,'forcedcalls',strcat(thissample(alpha,:),'myforcedscores.csv'))));
    microgliapositions=microgliapositions(:,11:19);
    allcellpositions=table2array(readtable(fullfile(mainfolderpath,'filteredallcellscounts',strcat(mynamethissample(alpha,:),'allcellsbarcodes.csv'))));
    allcellpositions=allcellpositions(:,3);
    maintypes=tdfread(fullfile(mainfolderpath,'kwanhocalls',strcat(thissample(alpha,:),'_allCells_major+PNsubtype.tsv')));
% DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
zeroImage=zeros(this.Height, this.Width);
%%
disp('starting plots')
    savepath=fullfile(mypath, 'analysis/newmosaics/');
    DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
    zeroImage=zeros(size(DAPI));dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
    homeo1neighbors=[];
    homeo2neighbors=[];
    immuneneighbors=[];
    inflammatoryneighbor=[];
    apeoneighbors=[];
    Ccr1neighbors=[];
    prolneighbors=[];
    allcells=maintypes.Oligo;
    allcells(:,2)=maintypes.CThPN;
    allcells(:,3)=maintypes.L6b_Subplate;
    allcells(:,4)=maintypes.L5_NP;
    allcells(:,5)=maintypes.INs;
    allcells(:,6)=maintypes.Astrocyte;  
    allcells(:,7)=maintypes.DL_CPN;
    allcells(:,8)=maintypes.L5_PT;
    allcells(:,9)=maintypes.Endothelia;  
    allcells(:,10)=maintypes.L5_CPN0x2FCStrPN;
    allcells(:,11)=maintypes.L20x2D3_CPN;
    allcells(:,12)=maintypes.L4_Stellate;    
    
    for k=1:7
        myImage=zeroImage;
        summed=reduceImage(zeroImage);
%         for i=1:size(modscore,1)
            nonzerovalues=fix(microgliapositions(:,9)).*(microgliapositions(:,k)>0);
            myImage(nonzerovalues(nonzerovalues>0))=1;
            nonzerovalues=nonzerovalues(find(nonzerovalues));
            neighbors=zeros(size(nonzerovalues,1),12);
            cellnumber=size(maintypes.Oligo,1);
            for thisMG=1:size(nonzerovalues,1)
                [row,column]=ind2sub(size(myImage),nonzerovalues(thisMG));
                distances=zeros(cellnumber,1);
                for thiscell=1:cellnumber
                    [cellrow,cellcolumn]=ind2sub(size(myImage),allcellpositions(thiscell,1));
                    distances(thiscell,1)=sqrt((cellrow-row)*(cellrow-row)+(cellcolumn-column)*(cellcolumn-column));
                end
                distances(distances>230)=0;
                if sum(distances)>0
                celltokeep=zeros(size(distances,1),1);
                celltokeep(find(distances))=1;
                cellIDs=allcells;
                mycellIDs=zeros(size(cellIDs));
                mycellIDs(find(cellIDs))=1;
                neighbors(thisMG,:)=sum(mycellIDs(logical(celltokeep),:),1);
                end
            end
    if k==1
            thisname='homeostatic1';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==2
            thisname='homeostatic2';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==3
            thisname='innate_immune';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==4
            thisname='inflammatory';
            writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==5
            thisname='ApoeHigh';
            ApoeHigh=reduceImage(imdilate(myImage,strel('disk', 30)));
           writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==6
        thisname='Ccr1High';
       writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));
        end
        if k==7
            thisname='proliferative';
        writematrix(neighbors,fullfile('/broad/hptmp/lbinan/jeffmicroglia/redo222/neighborhoods/',strcat(thissample(alpha,:), thisname,'.csv')));

        end
%         if k==1
%             thisname='homeostatic1';
%             homeo1=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+homeo1;
%             imwrite(homeo1,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         end
%         if k==2
%             thisname='homeostatic2';
%             homeo2=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+homeo2;
%             imwrite(homeo2,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         end
%         if k==3
%             thisname='innate_immune';
%             innate_immune=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+innate_immune;
%             imwrite(innate_immune,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         end
%         if k==4
%             thisname='inflammatory';
%             inflammatory=reduceImage(imdilate(myImage,strel('disk', 30)));
%             summed=summed+inflammatory;
%             imwrite(inflammatory,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         end
%         if k==5
%             thisname='ApoeHigh';
%             ApoeHigh=reduceImage(imdilate(myImage,strel('disk', 30)));
%             imwrite(ApoeHigh,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         end
%         if k==6
%         thisname='Ccr1High';
%         Ccr1High=reduceImage(imdilate(myImage,strel('disk', 30)));
%         summed=summed+Ccr1High;
%         imwrite(Ccr1High,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
% 
%         
%         end
%         if k==7
%             thisname='proliferative';
%         proliferative=reduceImage(imdilate(myImage,strel('disk', 30)));
%          summed=summed+proliferative;
%          imwrite(proliferative,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/prettymapsNEW/',strcat(thissample(alpha,:), thisname,'.png')));
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
%         imwrite(thisImage,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/kwanhoimages/',strcat(thissample(alpha,:), 'hex.png')));


end
% save('poolzscores.mat','poolzscores');