% clear;
clc;
close all;

dz=0.015;%���ÿռ䲽��
% dz=1;

Nx = 31;
Ny = 21;
ezx= (0:dz:(Nx-1)*dz);%���ÿռ������λ��
ezy= (0:dz:(Ny-1)*dz);%���ÿռ������λ��

for i=1:100%i��ȡֵ��ΧҪ����Ex��Hy�ľ���������Χ֮�ڣ�����ʱ��ȡ����
    Z=Ez(1+(i-1)*Ny:i*Ny,1:Nx);
%     Z=Hx(1+(i-1)*Ny:i*Ny,1:Nx);
    
    mesh(ezx,ezy,Z)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-64,28)%view E
%     view(0,0)
    %axis([0 (N-1)*dz 0 (N-1)*dz 0 0.1])%���������᷶Χ
    axis equal;
    axis([0 (Nx-1)*dz 0 (Ny-1)*dz 0 1])%���������᷶Χ
    
    pause(0.01)
end