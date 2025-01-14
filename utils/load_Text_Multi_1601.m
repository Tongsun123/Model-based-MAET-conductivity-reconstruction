function [Velocity, MAE, True_conductivity] = load_Text_Multi_1601(Velocity_filename, MAE_filename, Conductivity_filename)
% ���ܣ��������ݣ�ע�Ȿ����ֻ������������Ϊ1601������
% ���룺
% Code_path - �����·��
% Data_path - ���ݵ�·��
% Velocity_filename - �����ļ�������
% MAE_filename - �������ļ�����
% Conductivity_filename - �絼���ļ�����
% �����
% Velocity - �ٶ�
% MAE - MAE�ź�
% True_conductivity - ��ʵ�ĵ絼��
% Author��Sun Tong
% Data��2022-9-23
%%
% �����ٶ�����
Data = importdata(Velocity_filename);
data = Data.data;
Velocity = data(4802:end);
Velocity = circshift(-Velocity,-845);
% �������������
Data = importdata(MAE_filename);
data = Data.data;
t = data(1:1600);
MAE = data(4802:end);
% ����絼������
Data = importdata(Conductivity_filename);
data = Data.data;
True_conductivity = data(4802:end);
True_conductivity = True_conductivity / max(True_conductivity);