% clear;
clc;
close all;

dz=0.015;%设置空间步长

N = 21;
enz= (0:dz:(N-1)*dz);%设置空间采样点位置
eny= (0:dz:(N-1)*dz);%设置空间采样点位置

for i=1:200%i的取值范围要求在Ex和Hy的矩阵行数范围之内，代表时间取样点
    Z=Ez(1+(i-1)*N:i*N,1:N);
%     src(i)=Z(11,11);
    mesh(enz,eny,Z)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-64,28)%view E
    %axis([0 (N-1)*dz 0 (N-1)*dz 0 0.1])%设置坐标轴范围
    axis([0 (N-1)*dz 0 (N-1)*dz 0 1])%设置坐标轴范围
    pause(0.01)
end