%%
sound_velocity_curve_stand = ones(2000,1) * 1450;
distance_stand = cumsum(sound_velocity_curve_stand);
figure
plot(distance_stand)
% 曲线1对应有洞的
sound_velocity_curve1 = ones(2000,1);
sound_velocity_curve1(1:757) = sound_velocity_curve1(1:757) * 1450;
sound_velocity_curve1(758:869) = sound_velocity_curve1(758:869) * 1540;
sound_velocity_curve1(870:1119) = sound_velocity_curve1(870:1119) * 1450;
sound_velocity_curve1(1120:1299)= sound_velocity_curve1(1120:1299) * 1540;
sound_velocity_curve1(1300:end) = sound_velocity_curve1(1300:end) * 1450;
distance1 = cumsum(sound_velocity_curve1);
figure
plot(distance1)
% 曲线1对应没有洞的
sound_velocity_curve2 = ones(2000,1);
sound_velocity_curve2(1:756) = sound_velocity_curve2(1:756) * 1450;
sound_velocity_curve2(757:1271) = sound_velocity_curve2(757:1271) * 1540;
sound_velocity_curve2(1272:end) = sound_velocity_curve2(1272:end) * 1450;
distance2 = cumsum(sound_velocity_curve2);
figure
plot(distance2)
%
Threshold = 0.31;
slide_window = 60;
for j = 1:size(De_cond_hilbert,2)
    [peak_final, location_final] = Peak_find(De_cond_hilbert(:,j), Threshold, slide_window);
    number = size(peak_final,1);
    if number <= 3
        block = De_cond_hilbert(:,j);
        F = griddedInterpolant(distance2, block);                   % 形成grided算子F，我们可以对F采样得到均匀采集的k空间数据
        De_cond_hilbert(:,j) = F(distance_stand);
    elseif number >= 4
        block = De_cond_hilbert(:,j);
        F = griddedInterpolant(distance1, block);                   % 形成grided算子F，我们可以对F采样得到均匀采集的k空间数据
        De_cond_hilbert(:,j) = F(distance_stand);
    end
end

% block = De_cond_hilbert(:,check);
% F = griddedInterpolant(distance2, block);                   % 形成grided算子F，我们可以对F采样得到均匀采集的k空间数据
% block = F(distance_stand);
% figure
% plot(block,'b')
% hold on
% plot(De_cond_hilbert(:,check),'r')