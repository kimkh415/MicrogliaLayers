function tisdone=convertpixelIDbarcodes(mypath)
savepath=fullfile(mypath);
% load(fullfile(mypath,'positions2keep.mat'));
mkdir(savepath);
this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
zeroImage=zeros(this.Height, this.Width);
outputpath=fullfile(mypath);
barcodes=fix(table2array(readtable(fullfile(outputpath,'new_barcodes.csv'))));
essai=barcodes;
for i=1:size(essai,1)
    ind = sub2ind(size(zeroImage),essai(i,2),essai(i,1));
    essai(i,5)=ind;
end


% positions=ismember(essai(:,5),keeppositions);
% cleanedtable=essai(positions,:);
writematrix(essai,fullfile(savepath,'completedbarcodes.csv'));
% for i=1:7
%     thisZ=cleanedtable(cleanedtable(:,3)==i,:);
%     writematrix(thisZ,fullfile(savepath,strcat('barcodesforZ_',num2str(i),'.csv')));
% end
tisdone='done';
