
for i=1:10
    t(i)=995+i;
    newv(i) = new(t(i),t(i));
    oldv(i) = old(t(i),t(i));
end

plot(t,oldv,'r',t,newv,'b');
xlabel('n');
ylabel('�������');
legend('ԭ����','�ַ���','Location','northwest');