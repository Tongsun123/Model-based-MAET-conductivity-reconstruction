function Cond_model = Generate_conductivity_model(W, Nposition, Trans_band, Sample)
% ���ܣ��������й��ɴ��ĵ絼�ʺ���ģ��
% ���룺
% W - ��Ȩʸ��
% Nposition - �絼������λ��
% Trans_band - ���ɴ��ĵ���
% Sample - �ܲ��������Ŀ
% Authot��Sun Tong
% Date��2022-6-21
W = W(:);
Nposition = Nposition(:);
N = size(W,1);
Cond_model = zeros(Sample,1);                                 %�絼�ʵĳ���

for i=1:N
    step_function(1:Nposition(i)) = 0;
    step_function(Nposition(i)+1:Nposition(i)+Trans_band) = linspace(0,1,Trans_band); 
    step_function(end+1:Sample) = 1;
    step_function = step_function(:);
    Cond_model = Cond_model + W(i) * step_function;
end
end