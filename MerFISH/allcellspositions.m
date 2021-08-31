% thesepaths=['/broad/hptmp/lbinan/jeffmicroglia/Analysed/Run1/slice1side2/slice1side2/'];
upperCPN=[18,19,20,21,22,23];
L4stellate=[24,25,26];
deepcpn=[27,28,29];
scpn=[30,31,32,33,34,35,36,37];
chrpn=[38,39,40,41];
l6b=[42,43,44];
neurons=[45,46];
astrocytes=[47];
endothelium=[48,49,50];
oligo=[51:62];
interneuron=[63,64];
BAM=[65:67];
remember=[];
for alpha=1:size(thesepaths,1)
    mainfolderpath=['/broad/hptmp/lbinan/jeffmicroglia/'];
    mypath=thesepaths(alpha,:);
    allcellpositions=table2array(readtable(fullfile(mainfolderpath,'filteredallcellscounts',strcat(mypath(44:47),mypath(49:59),'allcellsbarcodes.csv'))));
    allcellpositions=allcellpositions(:,3);
    if beta~=5
    maintypes=tdfread(fullfile(mainfolderpath,'kwanhocalls',strcat(mypath(44:47),mypath(49:59),'_allCells_major+PNsubtype.tsv')));
    else
        maintypes=tdfread(fullfile(mainfolderpath,'kwanhocalls',strcat('Run3',mypath(49:59),'_allCells_major+PNsubtype.tsv')));
    end
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
    
    
    
    TupperCPN=[];
TL4stellate=[];
Tdeepcpn=[];
Tl5=[];
Tchrpn=[];
Tl6b=[];
TDL_CPN=[];
Tastrocytes=[];
Tendothelium=[];
TL202D3CPN=[];
Toligo=[];
Tinterneuron=[];
TL5CPN2FCStrPN=[];
TEndothelia=[];
TL5PT=[];
TIN=[];
TL4stellate=[];
    mypath=thesepaths(alpha,:);savepath=fullfile(mypath, 'analysis/newmosaics/');
    DAPI=imread(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
    zeroImage=zeros(size(DAPI));dapi=uint8(double(DAPI)/max(max(double(DAPI)))*255);
    mydapi=reduceImage(imbinarize(DAPI,'adaptive','sensitivity',0.3));
            if beta==1
                thissample=alpha;
            elseif beta==2
                thissample=alpha+2;
            else
                thissample=alpha+7;
            end
    for k=1:12
            myImage=zeroImage;
            nonzerovalues=fix(allcellpositions).*(allcells(:,k)>0);
            nonzerovalues=nonzerovalues(find(nonzerovalues));
            myImage(nonzerovalues(nonzerovalues>0))=1;
        if k==1
            thisname='oligo';
            Toligo=[Toligo;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==2
            thisname='Tchtpn';
            Tchrpn=[Tchrpn;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==3
            thisname='L6b_subplate';
            Tl6b=[Tl6b;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==4
            thisname='L5_NP';
            Tl5=[Tl5;scorepositions(thissample,zeroImage,myImage,mypath)];

        end
        if k==5
            thisname='INs';
            TIN=[TIN;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==6
            thisname='Astrocytes';
            Tastrocytes=[Tastrocytes;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==7
            IMthisname='DL_CPN';
            TDL_CPN=[TDL_CPN;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
        if k==8
            IMthisname='L5PT';
            TL5PT=[TL5PT;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
         if k==9
            IMthisname='Endothelia';
            TEndothelia=[TEndothelia;scorepositions(thissample,zeroImage,myImage,mypath)];
        end
         if k==10
            IMthisname='L5CPN2FCStrPN';
            TL5CPN2FCStrPN=[TL5CPN2FCStrPN;scorepositions(thissample,zeroImage,myImage,mypath)];
        end        
         if k==11
            IMthisname='L202D3CPN';
             TL202D3CPN=[TL202D3CPN;scorepositions(thissample,zeroImage,myImage,mypath)];
         end   
        if k==12
            IMthisname='L4stellate';
            TL4stellate=[TL4stellate;scorepositions(thissample,zeroImage,myImage,mypath)];
        end 

    end
    mysize=max([size(Tchrpn,1),size(Toligo,1),size(Tinterneuron,1),size(Tastrocytes,1)]);
    mydistances=zeros(mysize,7);
    mydistances(1:size(Toligo,1),1)=Toligo;
    mydistances(1:size(Tchrpn,1),2)=Tchrpn;
    mydistances(1:size(Tl6b,1),3)=Tl6b;
    mydistances(1:size(Tl5,1),4)=Tl5;
    mydistances(1:size(TIN,1),5)=TIN;
    mydistances(1:size(Tastrocytes,1),6)=Tastrocytes;
    mydistances(1:size(TDL_CPN,1),7)=TDL_CPN;   
    mydistances(1:size(TL5PT,1),8)=TL5PT;
    mydistances(1:size(TEndothelia,1),9)=TEndothelia;
    mydistances(1:size(TL5CPN2FCStrPN,1),10)=TL5CPN2FCStrPN;
    mydistances(1:size(TL202D3CPN,1),11)=TL202D3CPN;
    mydistances(1:size(TL4stellate,1),12)=TL4stellate;
    remember=[remember;mydistances];

end
writematrix(remember,fullfile('/broad/hptmp/lbinan/jeffmicroglia/Analysed/piaDist/',strcat('Run', num2str(beta),'Allcelldistances.csv')));

% homeo1=[1,0.8,0.4];
% homeo2=[0.8, 0, 0];
% apoe=[0.1216, 0.4706, 0.7059];
% innate_immune=[0.6510, 0.8078, 0.8902];
% ccr1=[0.4588, 0.4392,0.7020];
% inflammatory=[0.4, 0.651,0.1176];
% proliferaive=[0.9686,0.5059,0.7490];
