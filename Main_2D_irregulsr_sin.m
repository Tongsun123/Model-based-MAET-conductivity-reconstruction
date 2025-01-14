clear
clc
close all
RMSE = [];
SSIM_metric = [];
PSNR_metric = [];
for m = 1:1
%% Load filter parameter
tic
root_path = pwd;
addpath_config
tic
load 'High_pass_20MHz_0.30MHz.mat'
Sample = 1600;
Trans_band = 13;
Velocity_filename = '1MHz Velocity.txt';
MAE_filename = 'Model2_sin_1MHz_';
Conductivity_filename = '4_boundary_conductivity.txt';
Conductivity_filename2 = '6_boundary_conductivity.txt';
Conductivity_filename3 = '6_boundary_conductivity2.txt';
%% Load data and normalize
for i = -20:1:20
    MAE_filename_temp = [MAE_filename,num2str(i),'.txt'];
    [Velocity, MAE_temp, True_conductivity_interface4] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename);
    MAE(:,i+21) = MAE_temp;
end

Velocity = Velocity / max(Velocity);
[~, ~, True_conductivity_interface6] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename2);
[~, ~, True_conductivity_interface6_2] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename3);
MAE = MAE ./ max(MAE(:));
%% High-pass filter to remove low-frequency component
check = 21;
MAE_filtered = filtfilt(Num,1,MAE);
figure
plot(MAE_filtered(:,check))
%% Add the noise
SNR = 0.0;                                  % SNR Configuration 
bias = 0.0;                                 % Low-frequency Configuration
Noise = SNR * wgn(size(MAE_filtered,1),size(MAE_filtered,2),2)./max(wgn(size(MAE_filtered,1),size(MAE_filtered,2),2));
MAE_filtered = MAE_filtered + Noise * max(max(MAE_filtered)) + bias;                    % add noise
figure
plot(MAE_filtered(:,check))


figure
imagesc(1:41,14:99,MAE_filtered(400:1200,:))
colormap jet
axis off
axis equal
set(gca,'color','none');

alpha = 0.1;                % TGV alpha_1
beta = 0.1;                 % TGV alpha_0
nite = 30;                  % iteration
max_value = max(max(MAE));  % Normalization
MAE_filtered = TGV_denoise(MAE_filtered, alpha, beta, nite, max_value); % TGV denoising
%% Deconvolution
NSR = 1e3;
for i = 1:size(MAE_filtered,2)
    De_cond(:,i) = Wiener_Deconvolution_Multi(Velocity, MAE_filtered(:,i), NSR);
end

figure
plot(De_cond(:,check))

De_cond_hilbert = abs(hilbert(De_cond));

figure
plot(De_cond_hilbert(:,check))

De_cond_hilbert = De_cond_hilbert ./ max(max(De_cond_hilbert));

figure
imagesc(1:51,1:1601,De_cond_hilbert)
axis square

figure
imagesc(1:41,14:99,De_cond(400:1200,:))
colormap jet
axis off
axis equal
set(gca,'color','none');

[BW,maskedImage] = segmentImage_simulation1(De_cond_hilbert);
figure
imagesc(1:41,14:99,BW(400:1200,:))
axis off
axis equal
colormap gray

BW_grad = diff(BW);

figure
imagesc(1:41,14:99,BW_grad(400:1200,:))
axis off
axis equal
colormap gray
%% Locate the peak, establish a narrow threshold, and subsequently identify the authentic peak using a sliding window.
Threshold = 0.18;
slide_window = 80;
%
conductivity_matrix_C = [];
for j = 1:size(De_cond_hilbert,2)
    [peak_final, location_final] = Peak_find(De_cond_hilbert(:,j), Threshold, slide_window);
    %% Generate convolution matrix P
    Velocity_shift = circshift(-Velocity,-82);                  % 82
    convolution_matrix = zeros(size(Velocity_shift,1));
    for s = 1:size(convolution_matrix,1)
        convolution_matrix(s,:) = circshift(Velocity_shift,s);
    end
    signal = convolution_matrix * True_conductivity_interface4;
    %%
    step_function = zeros(size(convolution_matrix,1),1);
    conductivity_matrix = zeros(size(convolution_matrix,1),size(location_final,1));
    for i = 1:size(location_final,1)
        step_function(1:location_final(i)) = 0;
        step_function(location_final(i)+1:location_final(i)+Trans_band) = linspace(0,1,Trans_band); 
        step_function(location_final(i)+Trans_band+1:Sample) = 1;
        conductivity_matrix(:,i) = step_function;
    end
    
    conductivity_matrix_C = [conductivity_matrix_C,conductivity_matrix];
%% Conductivity reconstruction
    % Method 1
    lambda = 0;                                         % regularization parameter
    A = convolution_matrix * conductivity_matrix;       % Form A matrix (PC)
    
    AT = A';                                            %
    x_estimate = inv(AT * A+lambda * eye(size(AT,1))) * AT * MAE_filtered(:,j);     % Weighted computation by least squares
    x_estimate(end) = -sum(x_estimate(1:end-1));                                    % Apply constraints
    reconstruction_conductivity = conductivity_matrix * x_estimate;                 % Conductivity reconstruction
    reconstruction_conductivity_total(:,j) = reconstruction_conductivity;           % storage
    % Method 2
    x_estimate_polarity = x_estimate ./ abs(x_estimate);                            % Get polarity
    peak_correction = peak_final .* x_estimate_polarity;                            % Get the weight of deconvolution peak
    peak_correction(end) = -sum(peak_correction(1:end-1));                          % Apply constraints
    reconstruction_conductivity_correction = conductivity_matrix * peak_correction; % Corrected conductivity
    reconstruction_conductivity_correction_total(:,j) = reconstruction_conductivity_correction; % storage
%     figure
%     plot(reconstruction_conductivity,'b')
%     hold on
%     plot(True_conductivity_interface4,'r')

%     signal_recon = convolution_matrix * reconstruction_conductivity;
%     figure
%     plot(signal_recon,'b')
%     hold on
%     plot(MAE_filtered(:,j))
end
%% Visualization

figure
imshow(convolution_matrix,[])
colormap jet

figure
imagesc(1:2,14:99,conductivity_matrix_C(1:end,:))
axis off
axis square
colormap gray

AA = convolution_matrix * conductivity_matrix_C;
figure
imshow(AA(400:1200,:),[])
axis off
axis square
colormap jet


reconstruction_conductivity_total = reconstruction_conductivity_total / max(abs(reconstruction_conductivity_total(:)));
reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total / max(reconstruction_conductivity_correction_total(:));
figure
imagesc(1:41,14:99,reconstruction_conductivity_total(200:1400,:)) % 200:1400
colormap jet
axis off
axis equal

figure
imagesc(1:41,14:99,reconstruction_conductivity_correction_total(200:1400,:))
colormap jet
axis off
axis equal

reconstruction_conductivity_total_smooth = smoothn(reconstruction_conductivity_total,0.1,'robust');
figure
imagesc(1:41,14:99,reconstruction_conductivity_total_smooth(200:1400,:))
colormap jet
axis off
axis equal

reconstruction_conductivity_correction_total_smooth = smoothn(reconstruction_conductivity_correction_total,0.1,'robust');
figure
imagesc(1:41,14:99,reconstruction_conductivity_correction_total_smooth(200:1400,:))
colormap jet
axis off
axis equal

toc
MAE_check = convolution_matrix * reconstruction_conductivity_correction_total(:,check);
figure
plot(MAE_check / max(MAE_check),'b')
hold on
plot(MAE_filtered(:,check),'r')
%% Evaluation
Name = 'Model I_conductivity.txt';
[Temp] = importdata(Name);
z_direction = 1601;
y_direction = 41;
GT = zeros(z_direction,y_direction);
Temp_conductivity = Temp(:,3);
for i=1:z_direction
    for j = 1:y_direction
        GT(i,j) = Temp_conductivity(j+(i-1)*y_direction);
    end
end
GT = GT / max(GT(:));
figure
imagesc(1:41,14:99,GT(200:1400,:))
colormap jet
axis off
axis equal
mask = reconstruction_conductivity_correction_total_smooth > 0.1;
figure
imagesc(1:41,14:99,mask(200:1400,:))
colormap jet
axis off
axis equal
GT = GT .* mask;

RMSE_temp = norm(reconstruction_conductivity_correction_total_smooth - GT) / norm(GT);
SSIM_temp = ssim(reconstruction_conductivity_correction_total_smooth,GT);
PSNR_temp = psnr(reconstruction_conductivity_correction_total_smooth,GT);

RMSE = [RMSE,RMSE_temp];
SSIM_metric = [SSIM_metric,SSIM_temp];
PSNR_metric = [PSNR_metric,PSNR_temp];
% close all
end
toc