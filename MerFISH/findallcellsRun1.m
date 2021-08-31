function cellmask=findallcellsRun1(mypath,genelist)

this=imfinfo(fullfile(mypath,'GenerateMosaic','images','mosaic_DAPI_0.tif'));
zeroImage=zeros(this.Height, this.Width);
myImage=zeros(size(zeroImage));
load(fullfile(mypath,'positions2keep.mat'));


for i=1:size(genelist,2)
    for k=1:7
        myImage=myImage+generateImage(mypath,genelist(1,i),k);
    end
end
microglias=findcells(mypath);
myImage=myImage.*(ones(size(myImage))-microglias);
myCroppedImage=zeroImage;
myCroppedImage(keeppositions)=myImage(keeppositions);
myImage=imbinarize(myImage);
SE=strel('disk',3);
thisImage=imdilate(myImage,SE);
SE=strel('disk',3);
thisImage=imdilate(thisImage,SE);
SE=strel('disk',2);
thisImage=imerode(thisImage,SE);
SE=strel('disk',2);
thisImage=bwareaopen(imdilate(thisImage,SE),200);
thisImage=thisImage-bwareaopen(imdilate(thisImage,SE),1500);
SE=strel('disk',5);
thisImage=bwareaopen(imdilate(thisImage,SE),1000);
thisImage=thisImage-bwareaopen(imdilate(thisImage,SE),2500);
SE=strel('disk',6);
thisImage=imfill(bwareaopen(thisImage,800),'holes');imshow(thisImage)

cellmask=thisImage;
% rgb=cat(3,zeros(size(thisImage)),255*myImage,255*(imdilate(thisImage,SE)-thisImage));
% figure, imshow(rgb);