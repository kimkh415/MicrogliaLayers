function temp_barcodes=convertbarcodesCLUSTER(mypath)
%  mypath='/broad/hptmp/lbinan/jeffmicroglia/Analysed/Run2/slice3side2/slice3side2/';
outputpath=fullfile(mypath);
mkdir(outputpath);
% gene_list=[0;1;2;3];
barcodes=readtable(fullfile(mypath,'/barcodes.csv'));
conversion=readtable(fullfile(mypath,'micron_to_mosaic_pixel_transform.csv'));
temp_barcodes=barcodes(:,2:3);
temp_barcodes(:,3)=array2table(ones(size(temp_barcodes,1),1));
temp_barcodes=table2array(temp_barcodes);
newbarcodes=temp_barcodes*table2array(conversion).';
newbarcodes(:,3)=table2array(barcodes(:,4))+temp_barcodes(:,3);
newbarcodes(:,4)=table2array(barcodes(:,1))+temp_barcodes(:,3);
writematrix(newbarcodes,fullfile(outputpath,'new_barcodes.csv'));%x,y,z,gene
% myimage=zeros(size(mosaic_bit1_0));
