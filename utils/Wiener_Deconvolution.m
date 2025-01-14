function [Deconvoluted_Signal] = Wiener_Deconvolution(Velocity, MAEVol, NSR)
% Function£ºPerform deconvolution with Wiener filter
% Velocity - Velocity vector
% MAEVol - MAE signal vector
% NSR - regularizer
% Author£ºTong Sun
% Date£º2022-6-21
%% Normalization
Velocity = Velocity / max(Velocity);
MAEVol = MAEVol / max(MAEVol);
%% Transform time domain signal to frequency domain
Velocity_fft = fft(Velocity);
MAEVol_fft = fft(MAEVol);
%% Deconvolution
Deconvoluted_Signal_fft = MAEVol_fft.*conj(Velocity_fft)./(abs(Velocity_fft).^2+NSR);
Deconvoluted_Signal = ifft(Deconvoluted_Signal_fft);
%% Eliminate the zero-frequency offset after deconvolution.
Deconvoluted_Signal_Deviation = Deviation_Function(Deconvoluted_Signal(1:300));
Deconvoluted_Signal = Deconvoluted_Signal - Deconvoluted_Signal_Deviation;
Deconvoluted_Signal = Deconvoluted_Signal / max(abs(Deconvoluted_Signal));
end
