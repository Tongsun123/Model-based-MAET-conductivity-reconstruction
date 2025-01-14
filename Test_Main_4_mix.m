clear
clc
close all
%% Load filter parameter
root_path = pwd;
addpath_config
tic
load 'High_pass_40MHz_0.35MHz.mat'
Sample = 4001;
Trans_band = 20;
Velocity_filename = 'Velocity_4001.txt';
MAE_filename = 'MAE_4_boundary_mix.txt';
Conductivity_filename = '4_boundary_conductivity_mix.txt';
%% Load data and normalize
[Velocity, MAE, True_conductivity] = load_Text(Velocity_filename, MAE_filename, Conductivity_filename);
%% High-pass filter to remove low-frequency component
MAE_filtered = filtfilt(Num,1,MAE);
figure
plot(MAE_filtered)
%% Add the noise
SNR = 0.0;                                  % SNR Configuration 
Noise = SNR * wgn(size(MAE_filtered,1),size(MAE_filtered,2),2)./max(wgn(size(MAE_filtered,1),size(MAE_filtered,2),2));
MAE_filtered = MAE_filtered + Noise * max(max(MAE_filtered));                    % add noise
hold on 
plot(MAE_filtered)
hold off
%% Deconvolution
NSR = 5;
De_cond = Wiener_Deconvolution(Velocity, MAE_filtered, NSR);
De_cond = wiener2(De_cond);
figure
plot(De_cond)
De_cond_hilbert = abs(hilbert(De_cond));
figure
plot(De_cond_hilbert)
De_cond_hilbert = De_cond_hilbert ./ max(De_cond_hilbert);
%% Locate the peak, establish a narrow threshold, and subsequently identify the authentic peak using a sliding window.
Threshold = 0.2;
slide_window = 100;
[peak_final, location_final] = Peak_find(De_cond_hilbert, Threshold, slide_window);
%% Generate convolution matrix P
Velocity_shift = circshift(-Velocity,-360);
convolution_matrix = zeros(size(Velocity_shift,1));
for i = 1:size(convolution_matrix,1)
    convolution_matrix(i,:) = circshift(Velocity_shift,i);
end
signal = convolution_matrix * True_conductivity;
signal = signal / max(signal);
figure
plot(signal,'b')
hold on
plot(MAE,'r')
%%
step_function = zeros(size(convolution_matrix,1),1);
conductivity_matrix = zeros(size(convolution_matrix,1),size(location_final,1));
for i = 1:size(location_final,1)
    step_function(1:location_final(i)) = 0;
    step_function(location_final(i)+1:location_final(i)+Trans_band) = linspace(0,1,Trans_band); 
    step_function(location_final(i)+Trans_band+1:Sample) = 1;
    conductivity_matrix(:,i) = step_function;
end

% Verify whether the conductivity is generated correctly
% x = [0.5,1,-0.5,-1]';
% conductivity = conductivity_matrix * x;
% conductivity = conductivity / max(conductivity);
% figure
% plot(conductivity,'b')
% hold on
% plot(True_conductivity,'r')

% Verify whether the convolution matrix and conductivity matrix multiplied by 
% the correct x can obtain a signal consistent with the measured signal

% signal_2 = convolution_matrix * conductivity_matrix * x;
% signal_2 = signal_2 / max(signal_2);
% figure
% plot(signal_2,'b')
% figure
% plot(signal_2)
% hold on
% plot(MAE,'r')

%% Conductivity reconstruction
A = convolution_matrix * conductivity_matrix;
AT = A';
x_estimate = inv(AT * A) * AT * MAE_filtered;
x_estimate = x_estimate / max(x_estimate);
reconstruction_conductivity = conductivity_matrix * x_estimate;
reconstruction_conductivity = reconstruction_conductivity / max(reconstruction_conductivity);
figure
plot(reconstruction_conductivity,'b')
hold on
plot(True_conductivity,'r')
toc