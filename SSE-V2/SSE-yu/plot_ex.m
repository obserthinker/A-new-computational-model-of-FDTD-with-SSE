% clear;
clc;
close all;

dz=0.015;%���ÿռ䲽��

N = 21;
enz= (0:dz:(N-1)*dz);%���ÿռ������λ��
eny= (0:dz:(N-1)*dz);%���ÿռ������λ��

for i=1:200%i��ȡֵ��ΧҪ����Ex��Hy�ľ���������Χ֮�ڣ�����ʱ��ȡ����
    Z=Ez(1+(i-1)*N:i*N,1:N);
%     src(i)=Z(11,11);
    mesh(enz,eny,Z)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-64,28)%view E
    %axis([0 (N-1)*dz 0 (N-1)*dz 0 0.1])%���������᷶Χ
    axis([0 (N-1)*dz 0 (N-1)*dz 0 1])%���������᷶Χ
    pause(0.01)
end