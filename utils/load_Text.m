function [Velocity, MAE, True_conductivity] = load_Text(Velocity_filename, MAE_filename, Conductivity_filename)
% Function: Load data from txt file
% Input£º
% Velocity_filename - Name of velocity file
% MAE_filename - Name of MAE file
% Conductivity_filename - Name of Ground-truth conductivity

% Output£º
% Velocity
% MAE
% True_conductivity

% Load velocity datas
Data = importdata(Velocity_filename);
data = Data.data;
Velocity = data(12002:end);
Velocity = circshift(-Velocity,-1923);
Velocity = Velocity / max(Velocity);
% Load MAE signal
Data = importdata(MAE_filename);
data = Data.data;
t = data(1:4001);
MAE = data(12002:end);
MAE = MAE / max(MAE);
% Load Ground-truth conductivity
Data = importdata(Conductivity_filename);
data = Data.data;
True_conductivity = data(12006:end);
True_conductivity = True_conductivity / max(True_conductivity);