%%PLOT
clc; clear; clf('reset');
y = load('fort.21','-ascii');
plot(y(:,1),y(:,2));
hold on
for i=3:width(y)
    plot(y(:,1),y(:,i));
end
