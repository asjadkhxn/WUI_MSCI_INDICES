%%% horizon_average; used to construct independent variable, x
%%% Written by Yaojie Zhang  (��ҫ��); If you find any problem in this code, please send E-mail to yaojie.zhang@outlook.com. Thanks!
%%% First version: 2017/07/31; Last revised data:2017/7/31
function data_h=horizon_average_x(data,h)
%%% take average of data in the horizon of h.
%%data��ԭʼ������������,h��Ҫ�����ڣ�����h=5,��Ҫȡ������ľ�ֵ,
%%%   ��t�е�ȡ��ֵ������only use the data up to t. ǰh-1��Ϊnan

[n,m]=size(data); %n�죬m���ʲ�����

data_h=nan(n,m); %�ȹ��죬��Լ�㷨ʱ��

for j=1:m
    for i=h:n
        data_h(i,j)=mean(data(i-h+1:i,j));
    end    
end




