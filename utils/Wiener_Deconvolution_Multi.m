function [Deconvoluted_Signal] = Wiener_Deconvolution_Multi(Velocity, MAEVol, NSR)
% ���ܣ�ʹ��ά���˲���ʵʩ�����
% Velocity - ����ʸ��
% MAEVol - �������ź�ʸ��
% NSR - ������
% Author��Sun Tong
% Data��2022-6-21
%% ���ٺʹ������źŹ�һ��
Velocity = Velocity / max(Velocity);
%% ��ʱ���źű任��Ƶ��
Velocity_fft = fft(Velocity);
MAEVol_fft = fft(MAEVol);
%% ʵʩ�����
Deconvoluted_Signal_fft = MAEVol_fft.*conj(Velocity_fft)./(abs(Velocity_fft).^2+NSR);         %ά�����˲�
Deconvoluted_Signal = ifft(Deconvoluted_Signal_fft);                             %��fft
%% ��������������Ƶƫ��
Deconvoluted_Signal_Deviation = Deviation_Function(Deconvoluted_Signal(1:300));
Deconvoluted_Signal = Deconvoluted_Signal - Deconvoluted_Signal_Deviation;
end
