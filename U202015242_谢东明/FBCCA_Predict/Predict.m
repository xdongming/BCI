clc
clear
tic
freqs = 8:0.3:13.7;
phases = [];
for i = 1:5
    phases = [phases, (0:0.5:1.5)*pi];
end
N = 7;
a = 1.25;
b = 0.25;
wn = (1:N).^(-a) + b;
Nh = 4;
fs = 250;
num_temp = size(freqs, 2);
num_channel = 10;
num_trial = 22;
num_sample = 961;
Pred = zeros(num_trial, 10);
Template = zeros(num_temp, 2*Nh, num_sample);
names = {};
for i = 1:num_temp
    for n = 1:Nh
        t = 1:1:num_sample;
        Template(i,2*n-1,:) = sin(2*pi*n*freqs(i)/fs*t+phases(i));
        Template(i,2*n,:) = cos(2*pi*n*freqs(i)/fs*t+phases(i));
    end
end
for sub = 1:5
    for block = 1:2
        index = 2*(sub-1)+block;
        name = strcat('S',num2str(sub),'block',num2str(block));
        names{index} = name;
        load(strcat(name,'.mat'))
        start = find(data(11,:)==1);
        for j = 1:num_trial
            eeg = data(1:num_channel,start(j)+100:start(j)+4099);
            new_eeg = zeros(size(eeg,1),1000);
            for c = 1:size(eeg,1)
                new_eeg(c,:) = resample(eeg(c,:),fs,1000);
            end
            new_eeg = new_eeg(:, 0.16*fs:end);
            new_eeg = ideal_bandpassing(new_eeg, 2, 8, 49.8, fs);
            rho = zeros(num_temp,1);
            for k = 1:num_temp
                Rho_n = zeros(N,1);
                for n = 1:N
                    if n == 1
                        X = new_eeg';
                    else
                        X = ideal_bandpassing(new_eeg, 2, 8+4*(n-1), 49, fs)';
                    end
                    Y = squeeze(Template(k,:,:))';
                    [~, ~, r] = canoncorr(X, Y);
                    Rho_n(n) = r(1);
                end
                rho(k) = wn * (Rho_n.^2);
            end
            Pred(j,index) = find(rho == max(rho));
        end
    end
end
T = array2table(Pred, 'VariableNames', names);
filename = 'U202015242_Ð»¶«Ã÷.csv';
writetable(T, filename);