%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%程序功能：对走廊场景下的行人轨迹，结合PDR信息和TDOA信息，进行粒子滤波

%problem：
% (1)当步长偏移量deltaS不为0时，pdrTDOAMap会出现问题。
% (2)量测噪声设置得小，也会出现问题:可能是设置权值时，分母中会有量测噪声的缘故。

function [ errPdr,errPdrMap,errPdrTDOA,errPdrTDOAMap ] = PFTDOALoc()
clc;
close all;
clear;
% stepNum=100;
stepNum=36;
stepLength=zeros(1,stepNum) + 0.6; %每步步长为0.6米，步长真实值
theta=zeros(1,stepNum);%走直线，theta=0，真实值
theta(1,1:100) = 0;
% %设置声信号节点位置坐标
cornerLoc1 = [1.5;1.5];
cornerLoc2 = [22.5;1.5];
N = 200;    %粒子数量
% Q=0.1;
Q = 0.2;      %过程噪声（步长）。
Q1 = 0.01;       %过程噪声（角度）。
R = 0.3;      %量测噪声(TDOA测距误差)
staticStepLength = 0.8;     %初始步长设定值

T = length(stepLength) + 1; %状态数
X = zeros(2, T);    %存储系统状态，zeros零矩阵
TDOA = zeros(1, T);     %存储系统的观测量TDOA信息
pdr = zeros(2, T);     %PDR-Only
pdrMap = zeros(2,T);    %PDR+MAP
pdrTDOA = zeros(2,T);       %PDR+TDOA
pdrTDOAMap = zeros(2,T);    %PDR+TDOA+MAP

errPdr = zeros(1,T);    %纯PDR的误差
errPdrMap = zeros(1,T);     %PDR+MAP的误差
errPdrTDOA = zeros(1,T);    %PDR+TDOA的误差
errPdrTDOAMap = zeros(1,T);     %PDR+TDOA+MAP的误差

P = zeros(2, N);    %建立粒子群
w = zeros(N, 1);         %每个粒子的权重
pStepLength = zeros(1,N);   %每个粒子的步长
stepLengthEstimate = zeros(1,T);    %每一步步长的估计值

%%
%%初始化
startPointBiasX=0;
startPointBiasY=0;
X(:, 1) = [1.5+startPointBiasX; 1.5+startPointBiasY];     %初始系统状态,初始坐标
stepLengthEstimate(1) = staticStepLength;   %初始步长估计值

for i = 1 : N
%     P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2))+startPointBiasX;X(2,1)+wgn(1,1,10*log10(0.2))+startPointBiasY];
    P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2))-startPointBiasX;X(2,1)+wgn(1,1,10*log10(0.2))-startPointBiasY];
    pStepLength(i) = staticStepLength + wgn(1,1,10*log10(0.1));
%     dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离，向量的二范数：平方和开根号
%     w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重，高斯分布
    w(i) = 1;
end

pdr(:,1) = sum(P,2) / N;     %PDR的初始位置
pdrMap(:,1) = sum(P,2) / N;    
pdrTDOA(:,1) = sum(P,2) / N;       
pdrTDOAMap(:,1) = sum(P,2) / N;    

errPdr(1) = norm(X(:,1) - pdr(:,1));    %纯PDR的误差
errPdrMap(1) = norm(X(:,1) - pdrMap(:,1));     %PDR+MAP的误差
errPdrTDOA(1) = norm(X(:,1) - pdrTDOA(:,1));    %PDR+TDOA的误差
errPdrTDOAMap(1) = norm(X(:,1) - pdrTDOAMap(:,1));     %PDR+TDOA+MAP的误差

%%绘制可行区
map=[0,0;63,0;63,33;0,33];  %可行区，外边界
map1=[2,2;61,2;60,30;3,30]; %可行区，内边界
botSim1 = BotSim(map,[0,0,0]);
botSim = BotSim(map1,[0,0,0]);

%%
%%模拟运动轨迹，并设置量测值
for k = 2 : T
    X(:,k) = X(:,k-1) + stepLength(k-1) * [(cos(theta(k-1)));sin(theta(k-1))];
    TDOA(:,k) = norm(X(:,k) - cornerLoc1) - norm(X(:,k) - cornerLoc2) + wgn(1,1,10*log10(R));
end

%%
%%（1）PDR-Only、粒子滤波PDR+TDOA
for k = 2 : T   %运动时间为2~T时刻
    %预测
    theta(k-1) = theta(k-1) + wgn(1,1,10*log10(Q1))+0.1;             %考虑陀螺仪误差后的方向,发现方向偏差越小，则基于粒子滤波的动态步长估计算法，效果越好（因为在运动的垂直方向进行了约束）
    pdr(:, k) = pdr(:, k-1) + staticStepLength * [(cos(theta(k-1))); sin(theta(k-1))];  %使用固定步长情况下的PDR-Only 
%     stepLengthM = staticStepLength + wgn(1,1,10*log10(Q1));
    
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %每个粒子的步长
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        
        if mod(k,2) == 1    %k为奇数时
            dist = norm(P(:,i) - cornerLoc1) - norm(P(:,i) - cornerLoc2) - TDOA(k);
            w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);              
        end
    end
    
    %归一化权重
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %直接调用randomR函数（随机重采样算法）
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %另一种重采样方法
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则
        index = randi(N, 1);    %产生伪随机整数
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %得到新粒子
        pStepLength(i) = pStepLength(index);
    end
    pdrTDOA(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %计算误差
    errPdr(k) = norm(X(:,k) - pdr(:,k));
    errPdrTDOA(k) = norm(X(:, k) - pdrTDOA(:, k));     %粒子几何中心与系统真实状态的误差
end

%%（2）粒子滤波PDR+MAP
%重置
for i = 1 : N
    P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2));X(2,1)+wgn(1,1,10*log10(0.2))];
    pStepLength(i) = staticStepLength + wgn(1,1,10*log10(0.1));
    w(i) = 1;
end
in0 = zeros(1,N);
in1 = zeros(1,N);
in = zeros(1,N);
% preAdjacentPoint=zeros(2,T);
% distance = zeros(4,T);
stepLengthEstimate = zeros(1,T);    %每一步步长的估计值
stepLengthEstimate(1) = staticStepLength;

%开始滤波
for k = 2 : T
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %每个粒子的步长
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        
        %判断是否在可行区内，为1在可行区内，不为1则不在可行区内
        in0(i) = inpolygon(P(1,i),P(2,i),map(:,1),map(:,2)); %是否在外边界内
        in1(i) = inpolygon(P(1,i),P(2,i),map1(:,1),map1(:,2));%是否在内边界内
        in(i) = 0;
        if (in0(i) == 1) && (in1(i) ~= 1)
            in(i) = 1;
        else
            in(i) = 0;
        end        
    end
    
    %如果粒子在可行区，则将粒子权值设置为1，若粒子不在可行区，则将粒子权重设置为0
    for i = 1 : N
        if in(i) == 1
            w(i) = 1;           
        else
            w(i) = 0;
        end
    end    
    
    %归一化权重
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %直接调用randomR函数（随机重采样算法）
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %另一种重采样
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则
        index = randi(N, 1);    %产生伪随机整数
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %得到新粒子
        pStepLength(i) = pStepLength(index);
    end
    pdrMap(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %计算误差
    errPdrMap(k) = norm(X(:, k) - pdrMap(:, k));     %粒子几何中心与系统真实状态的误差    
end


%%（3）粒子滤波PDR+TDOA+MAP
%重置
for i = 1 : N
    P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2));X(2,1)+wgn(1,1,10*log10(0.2))];
    pStepLength(i) = staticStepLength + wgn(1,1,10*log10(0.1));
    w(i) = 1;
end
in0 = in0*0;
in1 = in1*0;
in = in*0;
% preAdjacentPoint=preAdjacentPoint*0;
% distance = distance*0;
stepLengthEstimate = stepLengthEstimate*0;    %每一步步长的估计值
stepLengthEstimate(1) = staticStepLength;

%开始滤波
for k = 2 : T
    %预测
%     theta(k-1) = theta(k-1) + wgn(1,1,10*log10(Q1*0.1));             %考虑陀螺仪误差后的方向,发现方向偏差越小，则基于粒子滤波的动态步长估计算法，效果越好（因为在运动的垂直方向进行了约束）
%     pdr(:, k) = pdr(:, k-1) + staticStepLength * [(cos(theta(k-1))); sin(theta(k-1))];  %使用固定步长情况下的PDR-Only 
    
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %每个粒子的步长
%         P(:, i) = P(:, i) + stepLengthM * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        PCenter(:,k) = sum(P,2)/N;   %当前粒子的几何中心
        %判断是否在可行区内，为1在可行区内，不为1则不在可行区内
        in0(i) = inpolygon(P(1,i),P(2,i),map(:,1),map(:,2)); %是否在外边界内
        in1(i) = inpolygon(P(1,i),P(2,i),map1(:,1),map1(:,2));%是否在内边界内
        in(i) = 0;
        if (in0(i) == 1) && (in1(i) ~= 1)
            in(i) = 1;
        else
            in(i) = 0;
        end        
    end
    
    %如果粒子在可行区，则将粒子权值设置为1，若粒子不在可行区，则将粒子权重设置为0
    for i = 1 : N
        if in(i) == 1
            if mod(k,2) == 1    %k为奇数时
                    dist = norm(P(:,i) - cornerLoc1) - norm(P(:,i) - cornerLoc2) - TDOA(k);
                    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);
            else
                 w(i) = 1;
            end            
        else
            w(i) = 0;
        end
    end    
    
    %归一化权重
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %直接调用randomR函数（随机重采样算法）
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %重采样（更新）?
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则
        index = randi(N, 1);    %产生伪随机整数
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %得到新粒子
        pStepLength(i) = pStepLength(index);
    end
    pdrTDOAMap(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %计算误差
    errPdrTDOAMap(k) = norm(X(:, k) - pdrTDOAMap(:, k));     %粒子几何中心与系统真实状态的误差    
end

% meanErrPdr = mean(errPdr);
% meanErrPdrMap = mean(errPdrMap);
% meanErrPdrTDOA = mean(errPdrTDOA);
% meanErrPdrTDOAMap = mean(errPdrTDOAMap);
end

% figure(2);
% plot(X(1,:),X(2,:),'r',pdr(1,:),pdr(2,:),'g',pdrMap(1,:), pdrMap(2,:),'b',pdrTDOA(1,:),pdrTDOA(2,:),'y',pdrTDOAMap(1,:),pdrTDOAMap(2,:),'lineWidth',1.2);
% % axis([0 200 0 40]);
% legend('True State', 'PDR-Only','PDR+MAP','PDR+TDOA','PDR+TDOA+MAP');
% hold on;
% botSim.drawMap();
% botSim1.drawMap();
% set(gca,'FontSize',12);
% xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
% 
% % 误差图
% figure(3);
% set(gca,'FontSize',12);
% plot(0:1:length(errPdr)-1,errPdr,'.-');
% hold on;
% plot(0:1:length(errPdrMap)-1,errPdrMap,'.-');
% hold on;
% plot(0:1:length(errPdrTDOA)-1,errPdrTDOA,'.-');
% hold on;
% plot(0:1:length(errPdrTDOAMap)-1,errPdrTDOAMap,'.-');
% xlabel('step', 'FontSize', 12);
% legend('pdr','pdrMap','pdrTDOA','pdrTDOAMap');
% title('The err');
% 
% % 概率密度分布图
% [f0,xi0] = ksdensity(errPdr(:));
% [f1,xi1] = ksdensity(errPdrMap(:));
% [f2,xi2] = ksdensity(errPdrTDOA(:));
% [f3,xi3] = ksdensity(errPdrTDOAMap(:));
% 
% figure(4);
% plot(xi0,f0);
% hold on;
% plot(xi1,f1);
% hold on;
% plot(xi2,f2);
% hold on;
% plot(xi3,f3);
% legend('pdr','pdrMap','pdrTDOA','pdrTDOAMap');
% title('PDF');
% 
% figure(5);
% cdfplot(abs(errPdr(:)));
% hold on;
% cdfplot(abs(errPdrMap(:)));
% hold on;
% cdfplot(abs(errPdrTDOA(:)));
% hold on;
% cdfplot(abs(errPdrTDOAMap(:)));
% legend('pdr','pdrMap','pdrTDOA','pdrTDOAMap');
% title('CDF');
%%
%随机重采样函数
function outIndex = randomR(weight)
L=length(weight);
outIndex=zeros(1,L);
u=unifrnd(0,1,1,L);
u=sort(u);
cdf=cumsum(weight);

i=1;
for j=1:L
    while((i <= L) && (u(i)<=cdf(j)))
        outIndex(i)=j;
        i=i+1;
    end
end
end
% end

