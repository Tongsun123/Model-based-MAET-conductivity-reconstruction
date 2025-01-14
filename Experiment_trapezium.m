clear
clc
close all
%%
root_path = pwd;
addpath_config
tic
% load 'High_pass_20MHz_0.30MHz.mat'
load 'Band_pass_20MHz_0.30-2MHz.mat'
Sample = 2000;
f_transducer = 1e6;
Trans_band = 13;
Velocity_filename = '1MHz Velocity.txt';
%% 导入数据并归一化
load 'Phantom_1.mat'
load 'H.mat'
time = y';
time = time(:,1);
MAE = x';
MAE(1:200,:) = MAE(1:200,:) * 0.01;
distance = time * 1450 * 1000;
selection = 415:1381;
Tra_position = 14:48;
Velocity = sin(2*pi*f_transducer*time);
figure
plot(time,Velocity)
rect_window = zeros(Sample,1);
rect_window(1:21) = 1;
Velocity = Velocity .* rect_window;
Velocity = circshift(Velocity,60);
Velocity = conv(Velocity,h);
Velocity = Velocity(1:Sample,1);
figure
plot(Velocity)
Velocity = circshift(Velocity,-300);
figure
plot(Velocity)
%% 高通滤波
check = 10;
% MAE_filtered = filtfilt(Num,1,MAE);
alpha = 0.5;
beta = 0.51;
nite = 10;
max_value = max(max(MAE));
MAE_filtered = TGV_denoise(MAE, alpha, beta, nite, max_value);
figure
plot(MAE_filtered(:,check))
%%
MAE_hilbert = abs(hilbert(MAE_filtered));
figure
plot(MAE_hilbert(:,check))
MAE_hilbert = MAE_hilbert ./ max(max(MAE_hilbert));
figure
imagesc(Tra_position,distance(selection),MAE_hilbert(selection,Tra_position))
axis off
axis equal
colormap jet
%% 找到峰值，并设置很小的阈值，再通过滑动窗找到真正的峰值
Threshold = 0.28;
slide_window = 60;
% 主循环
for j = 1:size(MAE_hilbert,2)
    [peak_final, location_final] = Peak_find(MAE_hilbert(:,j), Threshold, slide_window);
    %% 生成卷积矩阵
    Velocity_shift = circshift(Velocity,-165);          % 82
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
    %%
    step_function = zeros(size(convolution_matrix,1),1);
    conductivity_matrix = zeros(size(convolution_matrix,1),size(location_final,1));
    for i = 1:size(location_final,1)
        step_function(1:location_final(i)) = 0;
        step_function(location_final(i)+1:location_final(i)+Trans_band) = linspace(0,1,Trans_band); 
        step_function(location_final(i)+Trans_band+1:Sample) = 1;
        conductivity_matrix(:,i) = step_function;
    end
%% 重建电导率
    % 方法一
    lambda = 0;                                         % 正则化参数
    A = convolution_matrix * conductivity_matrix;       % 形成A矩阵(PC)
    AT = A';                                            % 得到A的转置
    x_estimate = inv(AT * A+lambda * eye(size(AT,1))) * AT * MAE_filtered(:,j);     % 最小二乘法得到权重
    x_estimate(end) = -sum(x_estimate(1:end-1));                                    % 施加约束
    reconstruction_conductivity = conductivity_matrix * x_estimate;                 % 重建电导率
    reconstruction_conductivity_total(:,j) = reconstruction_conductivity;           % 存储
    % 方法二
    x_estimate_polarity = x_estimate ./ abs(x_estimate);                            % 得到极性
    peak_correction = peak_final .* x_estimate_polarity;                            % 得到反卷积峰值的权重
    peak_correction(end) = -sum(peak_correction(1:end-1));                          % 施加约束
    reconstruction_conductivity_correction = conductivity_matrix * peak_correction; % 重建矫正后的电导率
    reconstruction_conductivity_correction_total(:,j) = reconstruction_conductivity_correction; % 存储
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
%% 可视化
figure
imshow(convolution_matrix,[])
colormap jet

reconstruction_conductivity_total = reconstruction_conductivity_total / max(abs(reconstruction_conductivity_total(:)));
reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total / max(reconstruction_conductivity_correction_total(:));
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_total(selection,Tra_position))
colormap jet
axis off
axis equal

figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_correction_total(selection,Tra_position))
colormap jet
axis off
axis equal

reconstruction_conductivity_total_smooth = smoothn(reconstruction_conductivity_total,0.1,'robust');
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_total_smooth(selection,Tra_position))
colormap jet
axis off
axis equal

reconstruction_conductivity_correction_total(:,14) = NaN;

reconstruction_conductivity_correction_total = reconstruction_conductivity_correction_total(selection,Tra_position);
reconstruction_conductivity_correction_total_smooth = smoothn(reconstruction_conductivity_correction_total,0.1,'robust');
figure
imagesc(Tra_position,distance(selection),reconstruction_conductivity_correction_total_smooth)
colormap jet
axis off
axis equal
toc