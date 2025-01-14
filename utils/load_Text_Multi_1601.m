function [Velocity, MAE, True_conductivity] = load_Text_Multi_1601(Velocity_filename, MAE_filename, Conductivity_filename)
% 功能：载入数据，注意本函数只能载入数据量为1601的数据
% 输入：
% Code_path - 代码的路径
% Data_path - 数据的路径
% Velocity_filename - 声速文件夹名字
% MAE_filename - 磁声电文件名字
% Conductivity_filename - 电导率文件名字
% 输出：
% Velocity - 速度
% MAE - MAE信号
% True_conductivity - 真实的电导率
% Author：Sun Tong
% Data：2022-9-23
%%
% 载入速度数据
Data = importdata(Velocity_filename);
data = Data.data;
Velocity = data(4802:end);
Velocity = circshift(-Velocity,-845);
% 载入磁声电数据
Data = importdata(MAE_filename);
data = Data.data;
t = data(1:1600);
MAE = data(4802:end);
% 载入电导率数据
Data = importdata(Conductivity_filename);
data = Data.data;
True_conductivity = data(4802:end);
True_conductivity = True_conductivity / max(True_conductivity);