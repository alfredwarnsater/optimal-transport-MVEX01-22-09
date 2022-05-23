clc;clf;clear

%f = @(x,y) x.^2 + (1+y).^2;
%fGrad = @(X) [2.*X(1); 2.*(1+X(2))];

f = @(x,y) x .* exp(-x.^2 - y.^2);
fGrad = @(X) [exp(-X(1).^2 - X(2).^2) .* (1 - 2.*X(1).^2); -2*exp(-X(1).^2 - X(2).^2).*X(1).*X(2)];

gridSize = 5;
resolution = 0.05;

[X,Y] = meshgrid(-gridSize:resolution:gridSize);
Z = f(X,Y);

%subplot(1,2,1)
%surf(X,Y,Z)
%grid on

%subplot(1,2,2)
contour(X,Y,Z,'LineWidth',2.5)
axis([-2 -0.05 -1.5 1.5])
title('Example of gradient descent', 'FontSize', 27)
xlabel('x', 'FontSize', 20)
ylabel('y', 'FontSize', 20)
grid on


%Gradient descent - - - - - - - - - - - - - - - - 
%start = [-2.2;1.2];
%maxItt = 100;
%h = 0.2;

start = [-1.5;0.88];
maxItt = 7;
h = 1;

[px, py] = gradient(Z);

cordValues = zeros(2, maxItt);
cordValues = [start,cordValues];

for i = 1:maxItt
    gradValues = -fGrad(cordValues(:,i));
    cordValues(:,i+1) = cordValues(:,i) + gradValues*h;
end

hold on
plot(cordValues(1,:),cordValues(2,:),'O-','color', 'r','LineWidth',2,'MarkerSize',12)
hold off


%hold on
%quiver(X,Y,px,py)
%hold off



