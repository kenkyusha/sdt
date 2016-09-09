close all
clear all
load('exampleData.mat');
n = 100;
x = linspace(-50,50,n);

y = [];
a = [];

w1 = [1 0]';
b1 = 0;

w2 = [0 1]';
b2 = 0;

ma = max(trainData);
mi = min(trainData);
n = size(trainData,1);

%normalization
 trainData = trainData - repmat(mi,n,1);
 trainData = trainData ./ repmat(ma,n,1);

b = SoftTree(trainData,trainTarget,trainData,trainTarget);
% b = SoftTree(trainData, trainTarget, validData, validTarget);
t = b.train();

x1 = linspace(min(trainData(:,1)),max(trainData(:,1)),n);
x2 = linspace(min(trainData(:,2)),max(trainData(:,2)),n);

p = [];
for k = 1:n
  for g = 1:n
   
    dataPoint = [x1(k) x2(g)];    
    p(k,g) = t.evaluate(dataPoint);
    
  end
end
figure
imagesc(x1,x2,p);axis xy;colorbar
% imagesc(x,x,a);axis xy;colorbar


hold on

plot(trainData(1:10,1),trainData(1:10,2),'ko');
hold on
plot(trainData(11:20,1),trainData(11:20,2),'ro');
legend;

% hold on
% 
% plot(trainData(1:400,1),trainData(1:400,2),'ko');
% hold on
% plot(trainData(401:800,1),trainData(401:800,2),'ro');




