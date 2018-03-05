function [C, imageBW_labeled] = labelBWImage(imageBW,filter_area)

fakeLabeledImage = bwlabel(imageBW, 8);     % Label each blob so we can make measurements of it
fakeMitoMeasurements = regionprops(fakeLabeledImage, imageBW, 'all');

allMitoAreas = [fakeMitoMeasurements.Area];
allowableAreaIndexes = allMitoAreas > filter_area; % Take the small objects.
keeperIndexes = find(allowableAreaIndexes);
keeperImage = ismember(fakeLabeledImage, keeperIndexes);
imageBW_labeled = bwlabel(keeperImage, 8);     % Label each blob so we can make measurements of it
C=bwconncomp(keeperImage,8);
end