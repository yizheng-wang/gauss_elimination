function [ x ] = gauss( A, b, order )
%�Զ���˹��ȥ�Լ������ֶ�ѡ����Ԫ�ĸ�˹��ȥ�����Ax=b
%    input�� A - ������
%            b - �������
%    output��order - ��Ԫ��ѡȡ������0���Զ���1������Ԫ��ȥ����-1����С����Ԫ������2��3��4...����Ԫ�Ĵ�С˳��
%                     �����ǵڶ��󣬵����󣬵��Ĵ�...��-2��-3��-4...��������Ԫ�Ĵ�С˳�������ǵڶ�С������С������С...
%Ĭ��order������0
if nargin < 3
    order = 0;
end
[rownum, colnum] = size(A);
Ab = [A, b]; % ƴ���������
% ��ȥ����
for j = 1 : colnum-1
            %ѡ��Ԫ
    if order ~= 0
            needorder = Ab(j:colnum, j);
            [~, I] = sort(needorder, 'descend');
    end
        
    if order >0
            getindex = I(order);
            Ab([getindex+j-1, j],:) = Ab([j, getindex+j-1],:);
    end
        
    if order <0
            getindex = I(end+order+1);
            Ab([getindex+j-1, j],:) = Ab([j, getindex+j-1],:);            
    end
    
    for i = j+1 : rownum
        l = Ab(i, j)/Ab(j, j);
        Ab(i, :) = Ab(i, :)- l * Ab(j, :);
    end
end
% �ش�����
x = zeros(colnum, 1);
for i = colnum : -1 : 1
    sum = 0 ;
    for m = i + 1 :colnum
        sum = sum + Ab(i, m) * x(m);
    end
    x(i) = (Ab(i, colnum + 1) - sum)/Ab(i, i);
end

end

