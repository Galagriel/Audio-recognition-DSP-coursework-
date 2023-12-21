function varargout = gui_audiorecognize(varargin)
% GUI_AUDIORECOGNIZE MATLAB code for gui_audiorecognize.fig
%      GUI_AUDIORECOGNIZE, by itself, creates a new GUI_AUDIORECOGNIZE or raises the existing
%      singleton*.
%
%      H = GUI_AUDIORECOGNIZE returns the handle to a new GUI_AUDIORECOGNIZE or the handle to
%      the existing singleton*.
%
%      GUI_AUDIORECOGNIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_AUDIORECOGNIZE.M with the given input arguments.
%
%      GUI_AUDIORECOGNIZE('Property','Value',...) creates a new GUI_AUDIORECOGNIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_audiorecognize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_audiorecognize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_audiorecognize

% Last Modified by GUIDE v2.5 07-Dec-2023 14:33:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_audiorecognize_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_audiorecognize_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_audiorecognize is made visible.
function gui_audiorecognize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_audiorecognize (see VARARGIN)

% Choose default command line output for gui_audiorecognize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_audiorecognize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_audiorecognize_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
duration = 2; % Recording duration (seconds)
recObj = audiorecorder(44100, 16, 2);
disp('正在录音……');
recordblocking(recObj, duration);
disp('录音已结束');
play(recObj); % Playback recorded audio
myRecording = getaudiodata(recObj);
fs =44100;
filename = 'D:\WYH\XJTU\作业\3-1数字信号处理\实验\temp\temp.wav'; % Temporary file to store recorded audio
audiowrite(filename,myRecording,44100);
filename = 'D:\WYH\XJTU\作业\3-1数字信号处理\实验\temp\temp.wav'; % Temporary file to store recorded audio
[y,fs]=audioread(filename);%使用 audioread 将数据读回 MATcLAB,fs为返回采样频率，y为读取的语音信号
T=1/fs; 
t=(0:length(y)-1)*T; %原始语音信号对应的采样时间序列
%以下进行归一化
y = y / max(abs(y));%幅度归一化到[-1,1]
%以下进行预加重y=x(i)-0.97*x(i-1)，预加重系数0.97
for i=2:88200
    x(i)=y(i)-0.95*y(i-1);
end
%以下进行双门限截断
Framelen=round(fs*0.02);%帧长为20ms
FrameLen=hann(Framelen);%选用海宁窗
inc=round(fs*0.01);%未重叠部分，即帧移为10ms
amp1 = 10;          %短时能量阈值
amp2 = 3;           %即设定能量的两个阈值。
zcr1 = 10;          %过零率阈值
zcr2 = 5;           %过零率的两个阈值，感觉第一个没有用到。
minsilence = 8;   %用无声的长度来判断语音是否结束
minlen  = 15;    %判断是语音的最小长度
status  = 0;      %记录语音段的状态
count   = 0;     %语音序列的长度
silence = 0;      %无声的长度
%计算过零率
tmp1  = enframe(x(1:end-1), Framelen,inc);%x为数据，Framelen为窗函数，inc是后一帧对前一帧的位移量
tmp2  = enframe(x(2:end), FrameLen,inc);%函数返回一个矩阵tmp,每行代表一个帧
signs = (tmp1.*tmp2)<0;
diffs = (tmp1 - tmp2)>0.005;
zcr   = sum(signs.*diffs,2);%对矩阵的行求和 %得到了信号各帧的过零率值，放到zcr矩阵中。
%计算短时能量
amp = sum((abs(enframe( x, FrameLen, inc))).^2, 2);%求出x各帧的能量值
%调整能量门限
amp1 = min(amp1, max(amp)/4);
amp2 = min(amp2, max(amp)/8);
%开始端点检测
for n=1:length(zcr)
    goto = 0;
    switch status
        case {0,1}                   % 0 = 静音, 1 = 可能开始
            if amp(n) > amp1          % 确信进入语音段
                x1 = max(n-count-1,1); % 记录语音段的起始点
                status  = 2;
                silence = 0;
                count   = count + 1;
            elseif amp(n) > amp2 || zcr(n) > zcr2 % 可能处于语音段
                status = 1;
                count  = count + 1;
            else                       % 静音状态
                status  = 0;
                count   = 0;
            end
        case 2                      % 2 = 语音段
            if amp(n) > amp2 ||zcr(n) > zcr2     % 保持在语音段
                count = count + 1;
            else                       % 语音将结束
                silence = silence+1;
                if silence < minsilence % 静音还不够长，尚未结束
                    count  = count + 1;
                elseif count < minlen   % 语音长度太短，认为是噪声
                    status  = 0;
                    silence = 0;
                    count   = 0;
                else                    % 语音结束
                    status  = 3;
                end
            end
        case 3
            break;
    end
end
count = count-silence/2;
x2 = x1 + count -1;              %记录语音段结束点
%后边的程序是找出语音端，然后用红线给标出来
%以下对归一化、预加重之后的原始信号进行截断
x=x(x1*inc:x2*inc);
filename = 'D:\WYH\XJTU\作业\3-1数字信号处理\实验\temp\temp1.wav'; % Temporary file to store recorded audio
audiowrite(filename,x,44100);
 %原始语音信号对应的采样时间序列
%set(handles.pushfbutton1, 'String', '开始录音'); % Update result display

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num_digits = 10; % 数字的数量（0~9）
mfcc_library = cell(10,1); % 创建存储MFCC特征的库，每个数字对应一个cell
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
                winlength=1024;
%                 win=hann(1);  
%                 %创建矩形窗
%                 win=rectwin(round(0.015 * fs));
%                 %创建布莱克曼窗
%                 win=blackman(round(0.015 * fs));
                coeffs = mfcc(y, fs,'NumCoeffs', 34,'Window',hamming(winlength),'OverlapLength', 0.5*winlength,'FFTLength',winlength);
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
for tens = 0:9
    digit_mfccs = mfcc_library{num2str(tens)}; % 获取当前数字的MFCC特征库 
    [test_y, test_fs] = audioread('D:\WYH\XJTU\作业\3-1数字信号处理\实验\temp\temp1.wav'); % Load recorded audio;
    test_coeffs = mfcc(test_y, test_fs, 'NumCoeffs', 34,'Window',hamming(winlength),'OverlapLength', 0.5*winlength,'FFTLength',winlength); % 提取测试语音的MFCC特征          
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
        distance = dtw(reshaped_coeffs, reshaped_digit_mfccs);
        if distance <= min_distance
            min_distance = distance;
            recognized_digit = digit; % 更新识别结果为当前数字
        end
    end
end
set(handles.pushbutton2, 'String', num2str(recognized_digit));

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pushbutton2, 'String', ''); % Clear the text in pushbutton2
