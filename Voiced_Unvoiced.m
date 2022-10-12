clear; close all;

run();

%Tìm ngưỡng bằng thuật toán nhị phân
function T = threshold(f, g)
    %Tìm Tmin Tmax từ f và g
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

    %Lặp nhị phân
    T = 1/2*(Tmin + Tmax);
    i = sum(f(f < T));
    p = sum(g(g > T));
    j = -1; q = -1;
    while i ~= j || p ~= q
        s = 1/length(f)*sum(f(f > T)) - 1/length(g)*sum(g(g < T));
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

%Tìm vùng chồng chéo (overlap area) truyền cho hàm threshold
function [T_ste, T_zcr] = getThreshold(timestamps, label, ste, zcr, frame_duration)
    f = []; g = [];     %ste uv và v
    f1 = []; g1 = [];   %zcr uv và v  
    f2 = [];            %ste sil
    j = 1;
    for i = 1:length(timestamps)
        if label{i} == "sil"
            while j*frame_duration <= timestamps(i) 
                f2(end+1) = ste(j);
                j = j + 1;
            end
        else
            if label{i} == "uv"    
                while j*frame_duration <= timestamps(i) 
                    f(end+1) = ste(j);
                    f1(end+1) = zcr(j);
                    j = j + 1;
                end
            else
                if label{i} == "v" 
                    while j*frame_duration <= timestamps(i)      
                        g(end+1) = ste(j);    
                        g1(end+1) = zcr(j);
                        j = j + 1;
                    end
                end
            end
        end
    end
    %Tìm ngưỡng
    T_ste = threshold(f, g);
    T_zcr = threshold(f1, g1);
end

%Tính STE của frame i
function result = STE(x, n, N)
    m = 0:N-1;
    result = sum(x(n-m).^2);
end

%Tính ZCR của frame i
function result = sgn(x)
    if x >= 0
        result = 1;
    else
        result = -1;
    end
end

function result = ZCR(x, n, N)
    result = 0;
    for m = 0:N-2
        result = result + abs(sgn(x(n-m))-sgn(x(n-m-1)));
    end
end

%Chuẩn hoá hàm thuộc tính về đoạn [-1;1]
function f = normalize(f, T, fmin, fmax)
    for i = 1:length(f)
        if f(i) >= T
            f(i) = (f(i)-T)/(fmax-T);
        else
            f(i) = (f(i)-T)/(T-fmin);
        end
    end
end

%Tính toán và vẽ đồ thị
function calAndPlot(wavFile, labFile, frame_duration)
    [x, Fs] = audioread(wavFile);

    %Số lượng mẫu trong 1 frame
    N = floor(frame_duration*Fs); % = frame_duration / T    
    frame_duration = N/Fs;

    %Tính số frame
    len = floor(length(x)/N);              
    N_rest = mod(length(x), N);      % N_rest: Số lượng mẫu của frame cuối
    if N_rest ~= 0
        len = len + 1;
    else
        N_rest = N;
    end

    %Tính ste và zcr
    ste = zeros(1,len); zcr = zeros(1, len);
    tmp = 0;
    for i = 1:len-1
        tmp = tmp + N;
        ste(i) = STE(x, tmp, N);
        zcr(i) = ZCR(x, tmp, N);
    end
    tmp = tmp + N_rest;
    ste(len) = STE(x, tmp, N_rest);
    zcr(len) = ZCR(x, tmp, N_rest);

    % Đọc dữ liệu các khung từ file
    timestamps = []; label = {};
    file = fopen(labFile);
    a = fscanf(file, "%*f%f");
    b = fscanf(file, "%s", 1);
    while b ~= "F0mean"
        timestamps(end+1) = a;
        label{end+1} = b;
        a = fscanf(file, "%*f%f");
        b = fscanf(file, "%s", 1);
    end
    fclose(file);

    [T_ste, T_zcr] = getThreshold(timestamps, label, ste, zcr, frame_duration);
    disp("Ngưỡng ste và zcr của audio " + wavFile + ": " + T_ste + ", " + T_zcr);
    
    %Chuẩn hoá ste và zcr về đoạn [-1;1] (khi đó ngưỡng sẽ đưa về 0)
    ste = normalize(ste, T_ste, min(ste), max(ste));
    zcr = normalize(zcr, T_zcr, min(zcr), max(zcr));

    %Tính VUijt
    VU = zeros(1, len);
    for i = 1:len
        if ste(i) - zcr(i) >= 0
            VU(i) = 1;
        end
    end
    
    %Tạo wave ste và zcr để plot
    ste_wave = zeros(1,length(x));
    zcr_wave = zeros(1,length(x));
    tmp = 0;
    for i = 1:len-1
        ste_wave(tmp + 1 : tmp + N) = ste(i);
        zcr_wave(tmp + 1 : tmp + N) = zcr(i);
        tmp = tmp + N;
    end
    ste_wave(tmp + 1 : tmp + N_rest) = ste(len);
    zcr_wave(tmp + 1 : tmp + N_rest) = zcr(len);

    %Vẽ đồ thị
    t = (1:length(x))/Fs;

    hold on;
    subplot(2, 1, 1); 
    plot(t, ste_wave,'r','LineWidth', 1); 
    hold on;
    plot(t, zcr_wave,'m','LineWidth', 1);
    plot(t, x, 'k');
    title('Trung gian'); xlabel('time(s)'); ylabel('x[n]');
    legend('Short-time Energy', 'Zero Cross Rate', 'Signal');
    hold off;

    subplot(2, 1, 2); 
    title('Cuối cùng'); xlabel('time(s)'); ylabel('x[n]');
    tmp = 0;
    for i = 1:len-1
        hold on;
        if VU(i) == 1
            p1 = plot(t(tmp + 1 : tmp + N), x(tmp + 1 : tmp + N),'b'); %voiced
        else
            p2 = plot(t(tmp + 1 : tmp + N), x(tmp + 1 : tmp + N),'c'); %unvoiced
        end
        tmp = tmp + N;
    end
    if VU(len) == 1
        plot(t(tmp + 1: tmp + N_rest), x(tmp + 1: tmp + N_rest),'b'); %voiced
    else
        plot(t(tmp + 1: tmp + N_rest), x(tmp + 1: tmp + N_rest),'c'); %unvoiced
    end

    border = [];    %Lưu vị trí biên
    %Vẽ biên phân biệt voiced-unvoiced
    for i = 1:len-1
        if VU(i) ~= VU(i+1)
            p3 = xline(i*N/Fs, 'g', 'LineWidth', 1);
            border(end+1) = i*N/Fs;
        end
    end  

    %Vẽ biên chuẩn phân biệt voiced-unvoiced
    for i = 1:length(timestamps)-1
        p4 = xline(timestamps(i), 'r','LineWidth', 1);
    end
    legend([p1 p2 p3 p4], {'Voiced', 'Unvoiced', 'My border', 'Standard border'});
    ylim([-1 1]);
    hold off;

    %disp(border); disp(timestamps);
    %Tính sai số RMSE giữa biên tìm được và biên chuẩn
    if length(border) > length(timestamps)
        minValue = length(timestamps);
    else
        minValue = length(border);
    end
    RMSE = 0;
    for i = 1:minValue
        RMSE = RMSE + (border(i)-timestamps(i))^2;
    end
    RMSE = sqrt(RMSE/minValue); 
    disp("Sai số RMSE giữa biên tìm được và biên chuẩn là: " + RMSE);
    disp("RMSE/SumTime = " + RMSE/t(end) + newline);
end

function run()
    frame_duration = 0.02;
    %Huấn luyện
%     figure('Name', 'Signal Phone FeMale');
%     calAndPlot("phone_F2.wav", "phone_F2.lab", frame_duration);
%     figure('Name', 'Signal Phone Male');
%     calAndPlot("phone_M2.wav", "phone_M2.lab", frame_duration);
%     figure('Name', 'Signal Studio Female');
%     calAndPlot("studio_F2.wav", "studio_F2.lab", frame_duration);
%     figure('Name', 'Signal Studio Male');
%     calAndPlot("studio_M2.wav", "studio_M2.lab", frame_duration);

    %Kiểm thử
    figure('Name', 'Signal Phone FeMale');
    calAndPlot("phone_F1.wav", "phone_F1.lab", frame_duration);
    figure('Name', 'Signal Phone Male');
    calAndPlot("phone_M1.wav", "phone_M1.lab", frame_duration);
    figure('Name', 'Signal Studio Female');
    calAndPlot("studio_F1.wav", "studio_F1.lab", frame_duration);
    figure('Name', 'Signal Studio Male');
    calAndPlot("studio_M1.wav", "studio_M1.lab", frame_duration);
end