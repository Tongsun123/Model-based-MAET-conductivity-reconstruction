clear
clc
close all

RMSE = [];
SSIM_metric = [];
PSNR_metric = [];
noise_mode = 'Lap';                             % Three kind of noise Gau, Lap, Uni
for m = 1:1
%% Load data and normalize
root_path = pwd;
addpath_config
tic
load 'High_pass_20MHz_0.30MHz.mat'
Sample = 1600;
Trans_band = 13;
Velocity_filename = '1MHz Velocity.txt';
MAE_filename = '';
Conductivity_filename = '4_boundary_conductivity.txt';
Conductivity_filename2 = '6_boundary_conductivity.txt';
Conductivity_filename3 = '6_boundary_conductivity2.txt';
%% Load data and normalize
for i = -25:1:25
    MAE_filename_temp = [MAE_filename,num2str(i),'.txt'];
    [Velocity, MAE_temp, True_conductivity_interface4] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename);
    MAE(:,i+26) = MAE_temp;
end
Velocity = Velocity / max(Velocity);
[~, ~, True_conductivity_interface6] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename2);
[~, ~, True_conductivity_interface6_2] = load_Text_Multi_1601(Velocity_filename, MAE_filename_temp, Conductivity_filename3);
MAE = MAE ./ max(MAE(:));
%% High-pass filter to remove low-frequency component
check = 35;
MAE_filtered = filtfilt(Num,1,MAE);
figure
plot(MAE_filtered(:,check))
%% Add the noise
SNR = 0.4;                                  % SNR Configuration 
bias = 0.0;                                 % …ËLow-frequency Configuration
switch(noise_mode)
    case 'Gau'
        Noise = SNR * wgn(size(MAE_filtered,1),size(MAE_filtered,2),2)./max(wgn(size(MAE_filtered,1),size(MAE_filtered,2),2));
        MAE_filtered = MAE_filtered + Noise * max(max(MAE_filtered)) + bias;                    % add noise
        figure
        hold on 
        plot(MAE_filtered(:,check))
        hold off
    case 'Lap'
        Noise_lap = SNR * ( 0 + 2.*randl(size(MAE_filtered,1),size(MAE_filtered,2)) ) ./ max(( 0 + 2.*randl(size(MAE_filtered,1),size(MAE_filtered,2)) ));
        MAE_filtered = MAE_filtered + Noise_lap * max(max(MAE_filtered)) + bias;
        figure
        hold on 
        plot(MAE_filtered(:,check))
        hold off
    case 'Uni'
        Noise_uni = SNR * (-1 + (1-(-1)).*rand(size(MAE_filtered,1),size(MAE_filtered,2))) ./ max(-1 + (1-(-1)).*rand(size(MAE_filtered,1),size(MAE_filtered,2)));
        MAE_filtered = MAE_filtered + Noise_uni * max(max(MAE_filtered)) + bias;
        figure
        hold on 
        plot(MAE_filtered(:,check))
        hold off
              
end
%%
alpha = 0.1;                % TGV alpha_1
beta = 0.1;                 % TGV alpha_0
nite = 30;                  % iteration
max_value = max(max(MAE));  % Normalization
MAE_filtered = TGV_denoise(MAE_filtered, alpha, beta, nite, max_value);
%% econvolution
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
%% Locate the peak, establish a narrow threshold, and subsequently identify the authentic peak using a sliding window.
Threshold = 0.25;
slide_window = 80;

for j = 1:size(De_cond_hilbert,2)
    [peak_final, location_final] = Peak_find(De_cond_hilbert(:,j), Threshold, slide_window);
    %% Generate convolution matrix P
    Velocity_shift = circshift(-Velocity,-82);          % 82
    convolution_matrix = zeros(size(Velocity_shift,1));
    for s = 1:size(convolution_matrix,1)
        convolution_matrix(s,:) = circshift(Velocity_shift,s);
    end
    
%     convolution_matrix = sparse(convolution_matrix);
%     [L,U] = ilu(convolution_matrix,struct('type','ilutp','droptol',1e-6));
%     [u,Flag,~,~,~] = gmres(convolution_matrix,MAE_filtered(:,20),[],1e-6,1000,L,U);
%     figure
%     plot(u)
%     figure
%     plot(MAE_filtered(:,20))
    
%     signal = convolution_matrix * True_conductivity_interface4;
    %%
    step_function = zeros(size(convolution_matrix,1),1);
    conductivity_matrix = zeros(size(convolution_matrix,1),size(location_final,1));
    for i = 1:size(location_final,1)
        step_function(1:location_final(i)) = 0;
        step_function(location_final(i)+1:location_final(i)+Trans_band) = linspace(0,1,Trans_band); 
        step_function(location_final(i)+Trans_band+1:Sample) = 1;
        conductivity_matrix(:,i) = step_function;
    end
%% Conductivity reconstruction
    % Method 1
    lambda = 0;                                         % regularization parameter
    A = convolution_matrix * conductivity_matrix;       % Form A matrix (PC)
    AT = A';
    x_estimate = inv(AT * A+lambda * eye(size(AT,1))) * AT * MAE_filtered(:,j); 
    x_estimate(end) = -sum(x_estimate(1:end-1));
    reconstruction_conductivity = conductivity_matrix * x_estimate;
    reconstruction_conductivity_total(:,j) = reconstruction_conductivity; 
    
    % Method 2
    x_estimate_polarity = x_estimate ./ abs(x_estimate);
    peak_correction = peak_final .* x_estimate_polarity;
    peak_correction(end) = -sum(peak_correction(1:end-1));
    reconstruction_conductivity_correction = conductivity_matrix * peak_correction;
    reconstruction_conductivity_correction_total(:,j) = reconstruction_conductivity_correction;
    
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

reconstruction_conductivity_total = reconstruction_conductivity_total / max(abs(reconstruction_conductivity_total(:)));
reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total / max(reconstruction_conductivity_correction_total(:));
figure
imagesc(1:51,14:99,reconstruction_conductivity_total(200:1400,:),[0,1])
colormap jet
axis off
axis equal

figure
imagesc(1:51,14:99,reconstruction_conductivity_correction_total(200:1400,:),[0,1])
colormap jet
axis off
axis equal

reconstruction_conductivity_total_smooth = smoothn(reconstruction_conductivity_total,0.1,'robust');
figure
imagesc(1:51,14:99,reconstruction_conductivity_total_smooth(200:1400,:),[0,1])
colormap jet
axis off
axis equal

% reconstruction_conductivity_correction_total(:,30) = NaN;
% reconstruction_conductivity_correction_total(:,44) = NaN;
% reconstruction_conductivity_correction_total(:,45) = NaN;

reconstruction_conductivity_correction_total_smooth = smoothn(reconstruction_conductivity_correction_total,0.1,'robust');
reconstruction_conductivity_correction_total_smooth = reconstruction_conductivity_correction_total_smooth / max(reconstruction_conductivity_correction_total_smooth(:));
figure
imagesc(1:51,14:99,reconstruction_conductivity_correction_total_smooth(200:1400,:),[0,1])
colormap jet
axis off
axis equal

MAE_check = convolution_matrix * reconstruction_conductivity_correction_total(:,check);
figure
plot(MAE_check / max(MAE_check),'b')
hold on
plot(MAE_filtered(:,check),'r')

%% RMSE
Name = 'Model II_conductivity.txt';
[Temp] = importdata(Name);
z_direction = 1601;
y_direction = 51;
GT = zeros(z_direction,y_direction);
Temp_conductivity = Temp(:,3);
for i=1:z_direction
    for j = 1:y_direction
        GT(i,j) = Temp_conductivity(j+(i-1)*y_direction);
    end
end
GT = GT / max(GT(:));
GT(701:750,18) = 0.6667;

figure
imagesc(1:51,14:99,GT(200:1400,:))
colormap jet
axis off
axis equal

reconstruction_conductivity_correction_total_smooth = circshift(reconstruction_conductivity_correction_total_smooth,[1,0]);
figure
imagesc(1:51,14:99,reconstruction_conductivity_correction_total_smooth(200:1400,:),[0,1])
colormap jet
axis off
axis equal

RMSE_temp = norm(reconstruction_conductivity_correction_total_smooth - GT) / norm(GT);
SSIM_temp = ssim(reconstruction_conductivity_correction_total_smooth,GT);
PSNR_temp = psnr(reconstruction_conductivity_correction_total_smooth,GT);

RMSE = [RMSE,RMSE_temp];
SSIM_metric = [SSIM_metric,SSIM_temp];
PSNR_metric = [PSNR_metric,PSNR_temp];
% close all
toc
end