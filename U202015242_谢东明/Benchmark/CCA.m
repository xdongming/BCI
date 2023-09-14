clc
clear
load('Freq_Phase.mat')
fs = 250;
subjects = 35;
num_trial = 40;
num_block = 6;
num_sample = 1216;
Pred = zeros(subjects, num_trial, num_block);
Nh = 5;
Template = zeros(num_trial, 2*Nh, num_sample);
for i = 1:num_trial
    for n = 1:Nh
        t = 1:1:num_sample;
        Template(i,2*n-1,:) = sin(2*pi*n*freqs(i)/fs*t+phases(i));
        Template(i,2*n,:) = cos(2*pi*n*freqs(i)/fs*t+phases(i));
    end
end
for num = 1:subjects
    name = strcat('S',num2str(num),'.mat');
    load(name)
    eeg=data([48 54:58 61:63], :, :, :);
    eeg_temp = eeg(:,0.64*fs :5.5*fs,:,:);
    [num_channel, num_sample, num_trial, num_block] = size(eeg_temp);
    new_eeg = zeros(num_channel, num_sample, num_trial, num_block);
    for i = 1:num_trial
        for j =1:num_channel
            for k =1:num_block
                new_eeg(j,:,i,k) = ideal_bandpassing(eeg_temp(j,:,i,k), 2, 8, 90, fs);
            end
        end
    end
    for j = 1:num_block
        for i = 1:num_trial
            rho = zeros(40,1);
            for k = 1:num_trial
                X = new_eeg(:,:,i,j)';
                Y = squeeze(Template(k,:,:))';
                [~, ~, r] = canoncorr(X, Y);
                rho(k) = r(1);
            end
            pred = find(rho == max(rho));
            if i == pred
                Pred(num, i, j) = 1;
            end
        end
    end
    fprintf('第%d代完成\n', num)
end
acc = sum(sum(sum(Pred)))./cumprod(size(Pred));
fprintf('准确率为%.2f%%\n',100*acc(end))