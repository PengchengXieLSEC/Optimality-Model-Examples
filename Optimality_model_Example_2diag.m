clear all
clc

% 可调变量
% f、Q
% x1：df、f、dQ、Q、x2、x2~、k => k0、k1
% w1：不等式条件第4项 <-> 反比例关系
% w2：k0_min、k1_max <-> 反比例关系，不等式条件第1、3、4项 <-> 正比例关系
% w3：x2、不等式条件第4项
% w4：x2~、k0_min、k1_max、不等式条件第2项
% eps：k0_min、k1_max、不等式条件第1、3项

% 调节步骤：
% 1.给定f、Q，预设一组w3、w4，调节x1，使计算得到的k中每个元素均小于-1（保证k0、k1都大于0）；
% 2.调节w2、eps，使k0、k1有取值范围，且不等式条件1、3项满足大小关系；
% 3.调节w3、w4，使不等条件第2项的值介于1、3项之间；如果调节后k0、k1的取值小于0，则返回第一步；
% 4.调节w1，使不等条件第4项的值大于第3项的值；至此，参数调节完毕。

% 参数配置
n = 2;

% x1 = 1;
x1 = [1;1.5];

w1 = 0.3;
w2 = 2;
w3 = 1;
w4 = 1;

% eps = 0.7;
eps = 0.8 / sqrt(n);
% eps = 1 / n;

% 计算f
fG = eye(n);
% fg = 1;
fg = [1;1];
fc = 0;
df = fG * x1 + fg;
f = 0.5 * x1' * fG * x1 + fg' * x1 + fc;

% 计算Q
QG = eye(n);
% Qg = 0;
Qg = [0;0];
Qc = 0;
dQ = QG * x1 + Qg;
Q = 0.5 * x1' * QG * x1 + Qg' * x1 + Qc;

% 计算x2和x2~
x2 = x1 - inv(fG + w3 * eye(n)) * df;
x2_t = x1 - inv(QG + w4 * eye(n)) * dQ;

% 计算k
k = (x2 - x1) ./ (x2_t - x2);

% 预算一下不等式条件(||·||F<=eps<1)是否存在解
A = inv(QG + w4 * eye(n)) * (eye(n) + diag(k));
B = inv(fG + w3 * eye(n)) * diag(k);
disp('B./A = ')
disp(B./A) % w2*A-w1*B < 0 => w2/w1 > B./A (1./A < 0)
X = eye(n) + w2 * inv(QG + w4 * eye(n)) * (eye(n) + diag(k)) - w1 * inv(fG + w3 * eye(n)) * diag(k);
disp(['系数矩阵的F范数为',num2str(norm(X,'fro'))])
disp(X)
disp(' ')

% 在如上条件下，考察k0和k1
k0_min = max(2 * (1 - eps) * (QG + w4 * eye(n)) / w2, [], "all");
k0_max = -1 * max(1 + k);% 按照定义，k0取此值
k1_min = -1 * min(1 + k);% 按照定义，k1取此值
k1_max = min(diag((QG + w4 * eye(n)) * (1 + eps) / w2), [], "all");
disp(['k0的范围：',num2str(k0_min),'<=k0<=',num2str(k0_max)])
disp(['k1的范围：',num2str(k1_min),'<=k1<=',num2str(k1_max)])
disp(' ')

if((k0_min > k0_max) || (k0_max < 0))
    disp('k0没有取值')
elseif ((k1_min > k1_max) || (k1_max < 0))
    disp('k1没有取值')
else
    % n = 1 与 n = 2, 对角
    disp('原条件中不等式从上到下四项依次为：')
    disp('第一项：')
    disp(k1_min * w2 / (1 + eps) * eye(n))
    disp('第二项：')
    disp(QG + w4 * eye(n))
    disp('第三项：')
    disp(k0_max * w2 / 2 / (1 - eps) * eye(n))
    disp('第四项：')
    disp(k0_max * w2 * (fG + w3 * eye(n)) / 2 / (1 + k1_min) / w1)

    % n = 2
%     disp('原条件中不等式从上到下四项依次为：')
%     disp('第一项：')
%     disp(k1_min * w2 / (1 + eps) * eye(n))
%     disp('第二项：')
%     disp(QG + w4 * eye(n))
%     disp('第三项：')
%     disp(k0_max * w2 / 2 / (1 - eps) * eye(n))
%     disp('第四项：')
%     disp(k0_max * w2 * (fG + w3 * eye(n)) / 2 / (1 + k1_min) / w1)
end
