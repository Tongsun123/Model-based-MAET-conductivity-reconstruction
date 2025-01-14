function [Deconvoluted_Signal] = Wiener_Deconvolution_Multi(Velocity, MAEVol, NSR)
% 功能：使用维纳滤波器实施反卷积
% Velocity - 声速矢量
% MAEVol - 磁声电信号矢量
% NSR - 正则化项
% Author：Sun Tong
% Data：2022-6-21
%% 声速和磁声电信号归一化
Velocity = Velocity / max(Velocity);
%% 将时域信号变换到频域
Velocity_fft = fft(Velocity);
MAEVol_fft = fft(MAEVol);
%% 实施反卷积
Deconvoluted_Signal_fft = MAEVol_fft.*conj(Velocity_fft)./(abs(Velocity_fft).^2+NSR);         %维纳逆滤波
Deconvoluted_Signal = ifft(Deconvoluted_Signal_fft);                             %逆fft
%% 消除反卷积后的零频偏移
Deconvoluted_Signal_Deviation = Deviation_Function(Deconvoluted_Signal(1:300));
Deconvoluted_Signal = Deconvoluted_Signal - Deconvoluted_Signal_Deviation;
end
