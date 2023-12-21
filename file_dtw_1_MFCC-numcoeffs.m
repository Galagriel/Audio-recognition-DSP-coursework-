% 创建一个存储不同numcoeffs下识别准确率的空数组
feature_counts = 12:40;
accuracies = zeros(size(feature_counts));
num_digits = 10; % 数字的数量（0~9）
mfcc_library = cell(10,1); % 创建存储MFCC特征的库，每个数字对应一个cell

% 循环不同的特征数量
for f = 1:numel(feature_counts)
    num_coeffs = feature_counts(f); % 当前的特征数量
    for tens = 0:9
        digit_mfccs = []; % 创建一个空的 MFCC 特征库用于存储当前数字的特征
    
        for ones = 0:20
            % 生成文件名，如果个位数小于10则在前面加0以匹配文件名格式
            if tens == 0
                file_name = sprintf('%03d.wav', ones);
            else
                file_name = sprintf('%d.wav', tens * 100+ ones);
            end
        
            % 构建完整文件路径
            file_path = fullfile("D:\WYH\XJTU\作业\3-1数字信号处理\实验\语音材料包\0-9\12月3号录音（换0）", file_name);
        
            % 检查文件是否存在，如果存在则读取该文件
            if exist(file_path, 'file') == 2
            
                [y, fs] = audioread(file_path);
                % 计算MFCC特征
                %创建汉宁窗
                winlength=512;
%                 disp(winlength);
% %                 win=hann();  
% %                 %创建矩形窗
% %                 win=rectwin(1024);
                coeffs = mfcc(y, fs,'NumCoeffs', num_coeffs,'Window',rectwin(winlength),'OverlapLength', 0.5*winlength,'FFTLength',winlength);
                % 将MFCC特征添加到当前数字的库中
                digit_mfccs = [digit_mfccs; coeffs];
            end
        end
    
        % 存储当前数字的MFCC特征库到mfcc_library中
        mfcc_library{num2str(tens)} = digit_mfccs;

        % 将当前数字的MFCC特征库保存为.mat文件到特定路径
        save_path = fullfile('D:\WYH\XJTU\作业\3-1数字信号处理\实验\result', ['mfcc_digit_' num2str(tens) '.mat']);
        save(save_path, 'digit_mfccs');
    end

    % 创建一个空的单元数组以存储所有测试结果
    results = cell(21*9, 4); % 21样本 * 9数字，每个样本包含“识别的数字”，“识别距离”，“文件名”，和“是否正确”

    index = 1; % 用于跟踪结果数组的索引
    correct_counts = zeros(1, 10); % 用于记录每个数字的正确识别数量
    total_counts = zeros(1, 10); % 用于记录每个数字的总样本数量
    
    for tens = 0:9
        digit_mfccs = mfcc_library{num2str(tens)}; % 获取当前数字的MFCC特征库 
        
        for ones = 0:20
            % 构建文件名和路径
            if tens == 0
                file_name = sprintf('%03d.wav', ones);
            else
                file_name = sprintf('%d.wav', tens * 100 + ones);
            end
            
            file_path = fullfile("D:\WYH\XJTU\作业\3-1数字信号处理\实验\testgroup\12月3号录音（换0）", file_name);
            
            if exist(file_path, 'file') == 2
                [test_y, test_fs] = audioread(file_path);
                test_coeffs = mfcc(test_y, test_fs, 'NumCoeffs', num_coeffs,'Window',rectwin(winlength),'OverlapLength', 0.5*winlength,'FFTLength',winlength); % 提取测试语音的MFCC特征
            
                [rows1, cols1, depth1] = size(test_coeffs);
                reshaped_coeffs = reshape(test_coeffs, rows1, cols1 * depth1);
            
                min_distance = Inf;
                recognized_digit = -1;
            
                for digit = 0:9
     
                    digit_mfccs = mfcc_library{num2str(digit)}; % 获取当前数字的MFCC特征库 
                    [rows2, cols2, depth2] = size(digit_mfccs);
                    reshaped_digit_mfccs = reshape(digit_mfccs, rows2, cols2 * depth2);   
                    [new_rows, new_cols] = size(reshaped_digit_mfccs);
                    if new_rows > rows1
                        reshaped_digit_mfccs = reshaped_digit_mfccs(1:rows1, :);
                    elseif new_rows < rows1
                        reshaped_digit_mfccs = [reshaped_digit_mfccs; zeros(rows1 - new_rows, new_cols)];
                    end
                
                    distance = dtw(reshaped_coeffs, reshaped_digit_mfccs,'absolute');
                %'squared'   'absolute'
                    if distance <= min_distance
                        min_distance = distance;
                        recognized_digit = digit; % 更新识别结果为当前数字
                    end
                end
            
                % 存储识别结果到结果数组
                results{index, 1} = recognized_digit;
                results{index, 2} = min_distance;
                results{index, 3} = file_name;
            
                % 计算每个数字的识别准确率
                if recognized_digit == tens
                    results{index, 4} = 1; % 标记为正确识别
                    correct_counts(tens + 1) = correct_counts(tens + 1) + 1;
                else
                    results{index, 4} = 0; % 标记为识别错误
                end
            
                total_counts(tens + 1) = total_counts(tens + 1) + 1;
            
                index = index + 1;
            end
        end
    end

    % 计算每个数字的识别准确率
    accuracy = correct_counts ./ total_counts;

    % 计算所有数字的平均识别准确率
    overall_accuracy = sum(correct_counts) / sum(total_counts);
    % 保存当前特征数量下的识别准确率
    accuracies(f) = overall_accuracy;
end
% 绘制折线图显示准确率随特征数量变化的趋势
plot(feature_counts, accuracies, 'o-');
xlabel('Numcoeffs');
ylabel('平均识别准确率(绝对距离)');
title('不同numcoeffs时0-9平均识别准确率(绝对距离)');
grid on;
% % 将结果存储到 Excel 文件
% result_file = "D:\WYH\XJTU\作业\3-1数字信号处理\实验\识别结果\频域.xlsx";
% xlswrite(result_file, results);
% disp('识别结果已保存到 Excel 文件中。');

% % 显示每个数字的识别准确率
% disp('每个数字的识别准确率：');
% for digit = 0:9
%     fprintf('数字 %d 的识别准确率为 %.2f%%\n', digit, accuracy(digit + 1) * 100);
% end
% 
% % 显示所有数字的平均识别准确率
% fprintf('所有数字的平均识别准确率为 %.2f%%\n', overall_accuracy * 100);


