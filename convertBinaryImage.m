function imageBW = convertBinaryImage(fname,frame)

image = imread(fname, frame);

%imageG = rgb2gray(image); don't need it if already grayscale
imageBW = im2bw(image,0.5);

% fill in holes in structures
imageBW = imfill(imageBW, 'holes');

end

