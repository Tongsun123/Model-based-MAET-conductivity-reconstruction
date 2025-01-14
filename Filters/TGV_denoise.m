function denoise = TGV_denoise(img, alpha, beta, nite, max_value)
% Function：Gaussian noise reduction of image using TGV
% 输入：
% img - Image to be denoiseds
%       alpha - alpha_1
%       beta - alpha_0
%       nite - Maximum number of iterations
%       max_value - Maximum values
% Ouput：
% denoise - Denoised image
img = img / max_value;
denoise = imtgvsmooth( img, alpha, beta, nite );
denoise = denoise * max_value;
end