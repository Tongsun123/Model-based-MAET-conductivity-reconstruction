function [BW,maskedImage] = segmentImage_pork_10_5(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 06-Oct-2024
%----------------------------------------------------


% Adjust data to span data range.
X = imadjust(X);

% 图割
foregroundInd = [12667 20667 21331 28667 29334 29338 29342 34667 35342 35346 42667 43350 43353 50667 50670 51353 51357 58670 59357 66670 67357 72670 73357 80670 88670 ];
backgroundInd = [21417 58746 58750 59399 72753 72757 72761 72840 72844 72848 72855 72882 72893 72908 80765 80867 80874 80919 88768 88772 88776 88942 96497 96780 96784 96787 96791 96957 102795 102802 102980 110806 110814 110821 110825 110829 110840 110848 110855 110859 110867 110870 110874 110878 111002 111433 119017 119033 119044 119055 127063 127070 127078 127085 127089 133093 133097 133101 133104 141116 141123 141127 141131 148942 149421 179372 ];
L = superpixels(X,1100);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

