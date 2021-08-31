% pathToYourImages='/this/is/your/path/';
 pathToYourImages='C:\Users\lbinan\Desktop\imagesfromjeff\neuronmaps';
imagesfiles=dir(fullfile(pathToYourImages,'*.png'));
for i=1:size(imagesfiles,1)
    myimage=imread(fullfile(imagesfiles(i).folder,imagesfiles(i).name));
    myname=imagesfiles(i).name;
    myname=myname(1:15);
    figure; imshow(2*myimage, []);
    roi = images.roi.Freehand;
    draw(roi);mask = createMask(roi);imshow(mask)
    stats=regionprops(mask,'PixelIdxList');
    L1=stats(1).PixelIdxList;
    name=convertCharsToStrings(fullfile(pathToYourImages,strcat(myname,'L1.mat')));
    save(name,'L1');
    figure; imshow(2*myimage, []);
    roi = images.roi.Freehand;
    draw(roi);mask = createMask(roi);imshow(mask)
    stats=regionprops(mask,'PixelIdxList');
    L2_4=stats(1).PixelIdxList;
    name=convertCharsToStrings(fullfile(pathToYourImages,strcat(myname,'L2_4.mat')));
    save(name,'L2_4');
    figure; imshow(2*myimage, []);
    roi = images.roi.Freehand;
    draw(roi);mask = createMask(roi);imshow(mask)
    stats=regionprops(mask,'PixelIdxList');
    L5=stats(1).PixelIdxList;
    name=convertCharsToStrings(fullfile(pathToYourImages,strcat(myname,'L5.mat')));
    save(name,'L5');
    figure; imshow(2*myimage, []);
    roi = images.roi.Freehand;
    draw(roi);mask = createMask(roi);imshow(mask)
    stats=regionprops(mask,'PixelIdxList');
    L6=stats(1).PixelIdxList;
    name=convertCharsToStrings(fullfile(pathToYourImages,strcat(myname,'L6.mat')));
    save(name,'L6');    
    close all
end