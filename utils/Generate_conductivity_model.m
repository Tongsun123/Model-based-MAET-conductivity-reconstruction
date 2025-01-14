function Cond_model = Generate_conductivity_model(W, Nposition, Trans_band, Sample)
% 功能：产生具有过渡带的电导率函数模型
% 输入：
% W - 加权矢量
% Nposition - 电导率跳变位置
% Trans_band - 过渡带的点数
% Sample - 总采样点的数目
% Authot：Sun Tong
% Date：2022-6-21
W = W(:);
Nposition = Nposition(:);
N = size(W,1);
Cond_model = zeros(Sample,1);                                 %电导率的长度

for i=1:N
    step_function(1:Nposition(i)) = 0;
    step_function(Nposition(i)+1:Nposition(i)+Trans_band) = linspace(0,1,Trans_band); 
    step_function(end+1:Sample) = 1;
    step_function = step_function(:);
    Cond_model = Cond_model + W(i) * step_function;
end
end