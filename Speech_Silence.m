clear; close all;

function T = threshold(f, g)
    if min(f) < min(g)
        Tmin = min(f);
    else
        Tmin = min(g);
    end

    if max(f) > max(g)
        Tmax = max(f);
    else
        Tmax = max(g);
    end

    T = 1/2*(Tmin + Tmax);
    i = sum(f(f < T));
    p = sum(g(g > T));
    j = -1; q = -1;
    while i ~= j || p ~= q
        s = 1/len(f)*sum(f(f > T)) - 1/len(g)*sum(g(g < T));
        if s > 0
            Tmin = T;
        else
            Tmax = T;
        end
        T = 1/2*(Tmin + Tmax);
        j = i; q = p;
        i = sum(f(f < T));
        p = sum(g(g > T));
    end
end
function main
    %[x, Fs] = audioread("phone_F2.wav");
    %for i = 1:

    [r, c] = size(frames);
    STE = zeros(r);
    for i = 1:r
        STE(i) = sum(pow(frames(i,:),2));
    end
    MA = sum(abs(x(n:1)));

    %Chuẩn hoá STE về đoạn [0:1]
    STE = (STE-min(STE))/(max(STE)-min(STE));


end
%Tìm ngưỡng bằng thống kê
function [T_ste, T_zcr] = getThreshold2(labFile, x, T, frame_duration)
    timstamp = [0.00,1.02,1.88,1.95,2.16,2.25,2.60,2.75,3.34,3.38,3.45,3.62,3.80,3.91,4.00,4.04,4.80];
    count_U = 0;
    count_V = 0;
    T_ste = 0; T_zcr = 0;
    mean_ste_U = 0;
    mean_ste_V = 0;
    mean_zcr_U = 0;
    mean_zcr_V = 0;

    for i = 1:length(timstamp)-1
        len_frame = (timstamp(i+1)-timstamp(i))/T;  %Số mẫu từ timstamp(i) - timstamp(i+1)
        n = timstamp(i)/T;                          %Số mẫu từ 0 - timstamp(i)
        if mod(i, 2) == 0
            mean_ste_V = mean_ste_V + STE(x, n, len_frame);
            mean_zcr_V = mean_zcr_V + ZCR(x, n, len_frame);
            count_V = count_V + 1;
        else
            mean_ste_U = mean_ste_U + STE(x, n, len_frame);
            mean_zcr_U = mean_zcr_U + ZCR(x, n, len_frame);
            count_U = count_U + 1;
        end
    end
    mean_ste_U = mean_ste_U/(count_U/frame_duration);
    mean_ste_V = mean_ste_V/(count_V/frame_duration);
    mean_zcr_U = mean_zcr_U/(count_U/frame_duration);
    mean_zcr_V = mean_zcr_V/(count_V/frame_duration);


end