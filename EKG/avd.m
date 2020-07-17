%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
clear
clc 
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://se.mathworks.com/help/signal/ug/find-periodicity-using-autocorrelation.html
%nahrání dat
%txt musí být ve stejné složce jako script
EKG = remove_noise(detrend(load ("mujSignal.txt")));
[EKG1_1,EKG1_2] = Graph_detrend("ecg_01.txt");
[EKG2_1,EKG2_2] = Graph_detrend("ecg_02.txt");
[EKG3_1,EKG3_2] = Graph_detrend("ecg_03.txt");
[EKG4_1,EKG4_2] = Graph_detrend("ecg_04.txt");
[EKG5_1,EKG5_2] = Graph_detrend("ecg_05.txt");
[EKG6_1,EKG6_2] = Graph_detrend("ecg_06.txt");
[EKG7_1,EKG7_2] = Graph_detrend("ecg_07.txt");
[EKG8_1,EKG8_2] = Graph_detrend("ecg_08.txt");
[EKG9_1,EKG9_2] = Graph_detrend("ecg_09.txt");
[EKG10_1,EKG10_2] = Graph_detrend("ecg_10.txt");

Corr_analysis(EKG,EKG10_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Funkce                               % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%funkce na vy?išt?ní signálu
function [ EKG1, EKG2 ] = Graph_detrend(file)
    EKG = load (file);
    %odstran?ní trendu grafu (normalizace hodnot)
    EKG(:,1) = detrend(EKG(:,1));
    EKG(:,2) = detrend(EKG(:,2));


    EKG1 = remove_noise(EKG(:,1));
    EKG2 = remove_noise(EKG(:,2));
end

% odstran?ní šumu ve funkci
function [removed_noise] = remove_noise(EKG)
    SP1_Y = EKG;
    hrubost = 2; % ..
    window = 7; 
    h = gausswin(hrubost*window+1)./window;
    Filt_EKG = zeros(size(SP1_Y));
    for i=1:length(SP1_Y)
        for j = -window:window
            if j>-i && j<(length(SP1_Y)-i+1)
                Filt_EKG(i) = Filt_EKG(i) + SP1_Y(i+j)*h(j+window+1);
            end
        end
    end

    removed_noise = Filt_EKG;
end

%funkce na korelaci a covaria?ní matici
function []= Corr_analysis(Template, Signal)

    %EKG = detrend(Template);
    Cov_Matrix = cov(Template,Signal);
    Correl = corrcoef(Template,Signal);
    fs = 200;   

    [autocor,lags] = xcorr(Signal,'coeff');
    [autocor2,lags2] = xcorr(Template,'coeff');
    
    figure('name','Auto-Correlation Graph - With all peaks');
    subplot(2,2,[1 3]); 
    plot(lags/fs,autocor,'b');
    grid on;
    hold on
    
    plot(lags2/fs,autocor2,'r');
    hold off;
    xlabel('Lag (Sec)');
    ylabel('Auto-Correlation');
    title('Auto-Correlation Graph - With all peaks');
    axis([-1 1 -0.1 1]);
    axis tight;
    
    h = subplot(2,2,2);
    title('Covarriance matrix');
    hPos = get(h, 'Position');
    t = uitable('Data', Cov_Matrix, 'ColumnName', {'Template', 'Signal'});
    t.Units = 'Normalized';
    t.Position = hPos;
    
    
    g = subplot(2,2,4);
    title('Correlation matrix');
    gPos = get(g, 'Position');
    ta = uitable('Data', Correl, 'ColumnName', {'Template', 'Signal'});
    ta.Units = 'Normalized';
    ta.Position = gPos;
    
    
end
