%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ܣ������ȳ����µ����˹켣�����PDR��Ϣ��TDOA��Ϣ�����������˲�

%problem��
% (1)������ƫ����deltaS��Ϊ0ʱ��pdrTDOAMap��������⡣
% (2)�����������õ�С��Ҳ���������:����������Ȩֵʱ����ĸ�л�������������Ե�ʡ�

function [ errPdr,errPdrMap,errPdrTDOA,errPdrTDOAMap ] = PFTDOALoc()
clc;
close all;
clear;
% stepNum=100;
stepNum=36;
stepLength=zeros(1,stepNum) + 0.6; %ÿ������Ϊ0.6�ף�������ʵֵ
theta=zeros(1,stepNum);%��ֱ�ߣ�theta=0����ʵֵ
theta(1,1:100) = 0;
% %�������źŽڵ�λ������
cornerLoc1 = [1.5;1.5];
cornerLoc2 = [22.5;1.5];
N = 200;    %��������
% Q=0.1;
Q = 0.2;      %������������������
Q1 = 0.01;       %�����������Ƕȣ���
R = 0.3;      %��������(TDOA������)
staticStepLength = 0.8;     %��ʼ�����趨ֵ

T = length(stepLength) + 1; %״̬��
X = zeros(2, T);    %�洢ϵͳ״̬��zeros�����
TDOA = zeros(1, T);     %�洢ϵͳ�Ĺ۲���TDOA��Ϣ
pdr = zeros(2, T);     %PDR-Only
pdrMap = zeros(2,T);    %PDR+MAP
pdrTDOA = zeros(2,T);       %PDR+TDOA
pdrTDOAMap = zeros(2,T);    %PDR+TDOA+MAP

errPdr = zeros(1,T);    %��PDR�����
errPdrMap = zeros(1,T);     %PDR+MAP�����
errPdrTDOA = zeros(1,T);    %PDR+TDOA�����
errPdrTDOAMap = zeros(1,T);     %PDR+TDOA+MAP�����

P = zeros(2, N);    %��������Ⱥ
w = zeros(N, 1);         %ÿ�����ӵ�Ȩ��
pStepLength = zeros(1,N);   %ÿ�����ӵĲ���
stepLengthEstimate = zeros(1,T);    %ÿһ�������Ĺ���ֵ

%%
%%��ʼ��
startPointBiasX=0;
startPointBiasY=0;
X(:, 1) = [1.5+startPointBiasX; 1.5+startPointBiasY];     %��ʼϵͳ״̬,��ʼ����
stepLengthEstimate(1) = staticStepLength;   %��ʼ��������ֵ

for i = 1 : N
%     P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2))+startPointBiasX;X(2,1)+wgn(1,1,10*log10(0.2))+startPointBiasY];
    P(:,i) = [X(1,1)+wgn(1,1,10*log10(0.2))-startPointBiasX;X(2,1)+wgn(1,1,10*log10(0.2))-startPointBiasY];
    pStepLength(i) = staticStepLength + wgn(1,1,10*log10(0.1));
%     dist = norm(P(:, i)-Z(:, 1));     %�����λ�����ľ��룬�����Ķ�������ƽ���Ϳ�����
%     w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ�أ���˹�ֲ�
    w(i) = 1;
end

pdr(:,1) = sum(P,2) / N;     %PDR�ĳ�ʼλ��
pdrMap(:,1) = sum(P,2) / N;    
pdrTDOA(:,1) = sum(P,2) / N;       
pdrTDOAMap(:,1) = sum(P,2) / N;    

errPdr(1) = norm(X(:,1) - pdr(:,1));    %��PDR�����
errPdrMap(1) = norm(X(:,1) - pdrMap(:,1));     %PDR+MAP�����
errPdrTDOA(1) = norm(X(:,1) - pdrTDOA(:,1));    %PDR+TDOA�����
errPdrTDOAMap(1) = norm(X(:,1) - pdrTDOAMap(:,1));     %PDR+TDOA+MAP�����

%%���ƿ�����
map=[0,0;63,0;63,33;0,33];  %����������߽�
map1=[2,2;61,2;60,30;3,30]; %���������ڱ߽�
botSim1 = BotSim(map,[0,0,0]);
botSim = BotSim(map1,[0,0,0]);

%%
%%ģ���˶��켣������������ֵ
for k = 2 : T
    X(:,k) = X(:,k-1) + stepLength(k-1) * [(cos(theta(k-1)));sin(theta(k-1))];
    TDOA(:,k) = norm(X(:,k) - cornerLoc1) - norm(X(:,k) - cornerLoc2) + wgn(1,1,10*log10(R));
end

%%
%%��1��PDR-Only�������˲�PDR+TDOA
for k = 2 : T   %�˶�ʱ��Ϊ2~Tʱ��
    %Ԥ��
    theta(k-1) = theta(k-1) + wgn(1,1,10*log10(Q1))+0.1;             %��������������ķ���,���ַ���ƫ��ԽС������������˲��Ķ�̬���������㷨��Ч��Խ�ã���Ϊ���˶��Ĵ�ֱ���������Լ����
    pdr(:, k) = pdr(:, k-1) + staticStepLength * [(cos(theta(k-1))); sin(theta(k-1))];  %ʹ�ù̶���������µ�PDR-Only 
%     stepLengthM = staticStepLength + wgn(1,1,10*log10(Q1));
    
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %ÿ�����ӵĲ���
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        
        if mod(k,2) == 1    %kΪ����ʱ
            dist = norm(P(:,i) - cornerLoc1) - norm(P(:,i) - cornerLoc2) - TDOA(k);
            w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);              
        end
    end
    
    %��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %ֱ�ӵ���randomR����������ز����㷨��
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %��һ���ز�������
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������
        index = randi(N, 1);    %����α�������
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�õ�������
        pStepLength(i) = pStepLength(index);
    end
    pdrTDOA(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %�������
    errPdr(k) = norm(X(:,k) - pdr(:,k));
    errPdrTDOA(k) = norm(X(:, k) - pdrTDOA(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����
end

%%��2�������˲�PDR+MAP
%����
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
stepLengthEstimate = zeros(1,T);    %ÿһ�������Ĺ���ֵ
stepLengthEstimate(1) = staticStepLength;

%��ʼ�˲�
for k = 2 : T
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %ÿ�����ӵĲ���
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        
        %�ж��Ƿ��ڿ������ڣ�Ϊ1�ڿ������ڣ���Ϊ1���ڿ�������
        in0(i) = inpolygon(P(1,i),P(2,i),map(:,1),map(:,2)); %�Ƿ�����߽���
        in1(i) = inpolygon(P(1,i),P(2,i),map1(:,1),map1(:,2));%�Ƿ����ڱ߽���
        in(i) = 0;
        if (in0(i) == 1) && (in1(i) ~= 1)
            in(i) = 1;
        else
            in(i) = 0;
        end        
    end
    
    %��������ڿ�������������Ȩֵ����Ϊ1�������Ӳ��ڿ�������������Ȩ������Ϊ0
    for i = 1 : N
        if in(i) == 1
            w(i) = 1;           
        else
            w(i) = 0;
        end
    end    
    
    %��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %ֱ�ӵ���randomR����������ز����㷨��
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %��һ���ز���
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������
        index = randi(N, 1);    %����α�������
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�õ�������
        pStepLength(i) = pStepLength(index);
    end
    pdrMap(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %�������
    errPdrMap(k) = norm(X(:, k) - pdrMap(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����    
end


%%��3�������˲�PDR+TDOA+MAP
%����
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
stepLengthEstimate = stepLengthEstimate*0;    %ÿһ�������Ĺ���ֵ
stepLengthEstimate(1) = staticStepLength;

%��ʼ�˲�
for k = 2 : T
    %Ԥ��
%     theta(k-1) = theta(k-1) + wgn(1,1,10*log10(Q1*0.1));             %��������������ķ���,���ַ���ƫ��ԽС������������˲��Ķ�̬���������㷨��Ч��Խ�ã���Ϊ���˶��Ĵ�ֱ���������Լ����
%     pdr(:, k) = pdr(:, k-1) + staticStepLength * [(cos(theta(k-1))); sin(theta(k-1))];  %ʹ�ù̶���������µ�PDR-Only 
    
    if k < 6
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,1:k-1));
    else
        tmpStepLengthEstimate = mean(stepLengthEstimate(:,k-5:k-1));
    end
    
    for i = 1 : N
        pStepLength(i) = tmpStepLengthEstimate + wgn(1,1,10*log10(0.01)); %ÿ�����ӵĲ���
%         P(:, i) = P(:, i) + stepLengthM * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        P(:, i) = P(:, i) + pStepLength(i) * [cos(theta(k-1)); sin(theta(k-1))] + wgn(2, 1, 10*log10(Q));
        PCenter(:,k) = sum(P,2)/N;   %��ǰ���ӵļ�������
        %�ж��Ƿ��ڿ������ڣ�Ϊ1�ڿ������ڣ���Ϊ1���ڿ�������
        in0(i) = inpolygon(P(1,i),P(2,i),map(:,1),map(:,2)); %�Ƿ�����߽���
        in1(i) = inpolygon(P(1,i),P(2,i),map1(:,1),map1(:,2));%�Ƿ����ڱ߽���
        in(i) = 0;
        if (in0(i) == 1) && (in1(i) ~= 1)
            in(i) = 1;
        else
            in(i) = 0;
        end        
    end
    
    %��������ڿ�������������Ȩֵ����Ϊ1�������Ӳ��ڿ�������������Ȩ������Ϊ0
    for i = 1 : N
        if in(i) == 1
            if mod(k,2) == 1    %kΪ����ʱ
                    dist = norm(P(:,i) - cornerLoc1) - norm(P(:,i) - cornerLoc2) - TDOA(k);
                    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);
            else
                 w(i) = 1;
            end            
        else
            w(i) = 0;
        end
    end    
    
    %��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
            w(i) = w(i) / wsum;
    end    
%     %ֱ�ӵ���randomR����������ز����㷨��
%     outIndex = randomR(w);
%     P(:,:) = P(:, outIndex);
    %�ز��������£�?
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������
        index = randi(N, 1);    %����α�������
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�õ�������
        pStepLength(i) = pStepLength(index);
    end
    pdrTDOAMap(:, k) = sum(P, 2) / N;
    stepLengthEstimate(k) = sum(abs(pStepLength)) / N;
    
    %�������
    errPdrTDOAMap(k) = norm(X(:, k) - pdrTDOAMap(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����    
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
% % ���ͼ
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
% % �����ܶȷֲ�ͼ
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
%����ز�������
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

