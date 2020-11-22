function [ x ] = gauss( A, b, order )
%自动高斯消去以及可以手动选择主元的高斯消去，求解Ax=b
%    input： A - 求解矩阵
%            b - 非齐次项
%    output：order - 主元的选取，其中0是自动，1是列主元消去法，-1是最小的主元，其余2，3，4...是主元的大小顺序
%                     依次是第二大，第三大，第四大...，-2，-3，-4...依次是主元的大小顺序，依次是第二小，第三小，第四小...
%默认order参数是0
if nargin < 3
    order = 0;
end
[rownum, colnum] = size(A);
Ab = [A, b]; % 拼成增广矩阵
% 消去过程
for j = 1 : colnum-1
            %选主元
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
% 回代过程
x = zeros(colnum, 1);
for i = colnum : -1 : 1
    sum = 0 ;
    for m = i + 1 :colnum
        sum = sum + Ab(i, m) * x(m);
    end
    x(i) = (Ab(i, colnum + 1) - sum)/Ab(i, i);
end

end

