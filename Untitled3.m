for pic = 60:118
    readpath = strcat('./result����pre/result/',num2str(pic),'_mask.png');
    img = imread(readpath);
    img_new = imresize(img,[1080 1920]);
    savepath = strcat('./result_new/',num2str(pic),'_mask.png');
    imwrite(img_new,savepath)
end