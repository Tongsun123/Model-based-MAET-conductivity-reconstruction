clear
clc
close all
%%
root_path = pwd;
addpath_config
tic
% load 'High_pass_20MHz_0.30MHz.mat'
load 'Band_pass_20MHz_0.15-1.6MHz.mat'
Sample = 2000;
f_transducer = 1e6;
Trans_band = 40; 
Velocity_filename = '1MHz Velocity.txt';
%%
load 'sin_3cycle_pork_800mV_2000ave_10_5.mat'
load 'H.mat'
time = x(1:2000,:);
time = time(:,1);
MAE = y(1:2000,:);
MAE(1:300,:) = MAE(1:300,:) * 0.01;
MAE(1500:end,:) = MAE(1500:end,:) * 0.1;
distance = time * 1450 * 1000;
selection = 415:1519;
Tra_position = 1:50;
Velocity = sin(2*pi*f_transducer*time);
figure
plot(time,Velocity)
rect_window = zeros(Sample,1);
rect_window(1:60) = 1;
Velocity = Velocity .* rect_window;
Velocity_H = circshift(Velocity,60);
Velocity_H = filtfilt(Num,1,Velocity_H);
Velocity_H = Velocity_H(1:Sample,1);
figure
plot(Velocity_H)
Velocity_H = circshift(-Velocity_H,-30);
figure
plot(Velocity_H)
figure
imagesc(Tra_position,distance(selection),MAE(selection,Tra_position))
axis off
axis equal
colormap jet
att_curve = exp(0.005*(distance-30));
figure
plot(att_curve)

MAE = MAE .* att_curve;
figure
imagesc(Tra_position,distance(selection),MAE(selection,Tra_position))
axis off
axis equal
colormap jet
%%
check = 40;
figure
plot(y(:,check))
figure
plot(MAE(:,check))

MAE(:,38:52) = MAE(:,38:52) * 1.4;
MAE(:,53:60) = MAE(:,53:60) * 1;

% MAE = filtfilt(Num,1,MAE);
figure 
plot(MAE(:,check))
figure
imagesc(Tra_position,distance(selection),MAE(selection,Tra_position))
axis off
axis equal
colormap jet

alpha = 0.4;
beta = 0.41;
nite = 50;
max_value = max(max(MAE));
MAE_filtered = TGV_denoise(MAE, alpha, beta, nite, max_value);
% MAE_filtered = MAE;
figure
plot(MAE_filtered(:,check))
MAE_filtered = MAE_filtered / max(MAE_filtered(:));
figure
imagesc(Tra_position,distance(selection),MAE_filtered(selection,Tra_position))
axis off
axis equal
colormap jet
%%
NSR = 5e3;
for i = 1:size(MAE_filtered,2)
    De_cond(:,i) = Wiener_Deconvolution_Multi(Velocity_H, MAE_filtered(:,i), NSR);
end
figure
plot(De_cond(:,check))
De_cond_hilbert = abs(hilbert(De_cond));
figure
plot(De_cond_hilbert(:,check))
De_cond_hilbert = De_cond_hilbert ./ max(max(De_cond_hilbert));

figure
imagesc(Tra_position,distance(selection),De_cond_hilbert(selection,Tra_position))
axis off
axis equal
colormap jet
%%
[BW,maskedImage] = segmentImage_pork_10_5(De_cond_hilbert);
figure
imagesc(Tra_position,distance(selection),BW(selection,Tra_position))
axis off
axis equal
colormap jet
%%
% sound_velocity_correction
% curve = ones(80,1)*1.2;
% curve(1:12) = 2;
% curve(21:27) = 1.5;
% curve(38:44) = 1.5;
% curve(67:74) = 2;
% for i = 1:size(MAE,2)
%     De_cond_hilbert(1250:1400,i) = De_cond_hilbert(1250:1400,i) .* curve(i);
% end
% figure
% imagesc(1:80,1:145,De_cond_hilbert(:,1:80))
% colormap jet
% axis square
%%
Threshold = 0.1;
slide_window = 200;

BW_grad = diff(BW);

figure
imagesc(Tra_position,distance(selection),BW_grad(selection,Tra_position))
axis off
axis equal
colormap jet

for j = 1:60
    j  
    %
    [peak_final_pre_section, location_final_pre_section] = Peak_find(De_cond_hilbert(:,j), Threshold, slide_window);
    
    
    rising = find(BW_grad(:,j)>0);
    falling = find(BW_grad(:,j)<0);
    location_final = zeros(size(rising));
    peak_final = zeros(size(rising));
    
    for m = 1:size(rising,1)
        for n = 1:size(location_final_pre_section,1)
            Flag = location_final_pre_section(n)>rising(m) && location_final_pre_section(n)<falling(m);
            if Flag == 1
                location_final(m) = location_final_pre_section(n);
                Flag = 0;
            end
        end

    end
    
    for i = 1:size(location_final,1)
        if location_final(i)==0
            location_final(i) = round((rising(i) + falling(i))/2);
        end
        peak_final(i) = De_cond_hilbert(location_final(i),j);    
    end
    
    Velocity_shift = circshift(-Velocity_H,-82);          % 82
    convolution_matrix = zeros(size(Velocity_shift,1));
    for s = 1:size(convolution_matrix,1)
        convolution_matrix(s,:) = circshift(Velocity_shift,s);
    end
    
    step_function = zeros(size(convolution_matrix,1),1);
    conductivity_matrix = zeros(size(convolution_matrix,1),size(location_final,1));
    for i = 1:size(location_final,1)
        step_function(1:location_final(i)) = 0;
        step_function(location_final(i)+1:location_final(i)+Trans_band) = linspace(0,1,Trans_band); 
        step_function(location_final(i)+Trans_band+1:Sample) = 1;
        conductivity_matrix(:,i) = step_function;
    end
     %
    lambda = 0;                                         %
    A = convolution_matrix * conductivity_matrix;       %
    AT = A';                                            %
    x_estimate = inv(AT * A+lambda * eye(size(AT,1))) * AT * MAE_filtered(:,j);     %
    if x_estimate(1)<0
        x_estimate = x_estimate * (-1);
    end
    x_estimate(end) = -sum(x_estimate(1:end-1));                                    %
    reconstruction_conductivity = conductivity_matrix * x_estimate;                 %
    reconstruction_conductivity_total(:,j) = reconstruction_conductivity;           %
    %
    x_estimate_polarity = x_estimate ./ abs(x_estimate);                            %
    peak_correction = peak_final .* x_estimate_polarity;                            %
    peak_correction(end) = -sum(peak_correction(1:end-1));                          %
    reconstruction_conductivity_correction = conductivity_matrix * peak_correction; %
    reconstruction_conductivity_correction_total(:,j) = reconstruction_conductivity_correction; %
end
%%
figure
imshow(convolution_matrix,[])
colormap jet

reconstruction_conductivity_total = reconstruction_conductivity_total / max(abs(reconstruction_conductivity_total(:)));
reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total / max(reconstruction_conductivity_correction_total(:));
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_total(selection,Tra_position))
axis off
axis equal
colormap jet

figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_correction_total(selection,Tra_position),[0,1])
axis off
axis equal
colormap jet

reconstruction_conductivity_correction_total(:,42:44) = NaN;

reconstruction_conductivity_total_smooth = smoothn(reconstruction_conductivity_total,0.1,'robust');
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_total_smooth(selection,Tra_position))
axis off
axis equal
colormap jet

reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total(selection,Tra_position);
reconstruction_conductivity_correction_total_smooth = smoothn(reconstruction_conductivity_correction_total,0.1,'robust');
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_correction_total_smooth,[0,1])
axis off
axis equal
colormap jet

reconstruction_conductivity_correction_total_smooth = reconstruction_conductivity_correction_total_smooth / max(reconstruction_conductivity_correction_total_smooth(:));
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_correction_total_smooth,[0,1])
axis off
axis equal
colormap jet

% MAE_check = convolution_matrix * reconstruction_conductivity_correction_total(:,check);
% figure
% plot(MAE_check / max(MAE_check),'b')
% hold on
% plot(MAE_filtered(:,check),'r')