k = 1:K_max;
sz = size(A_script);
N = sz(1);
sz1 = size(A);
n = sz1(1);
x_R = x_r;
for i = 2:N
    x_R = [x_R; x_r];
end

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
for i = 1:N-1
    x1 = [x1; x(1 + n*i,:)];
    x2 = [x2; x(2 + n*i, :)];
    x3 = [x3; x(3 + n*i, :)];
end


subplot(4,1,1);
plot(k,x1, 'LineWidth', 1)
hold on
plot(k,x_r(1,:), 'LineWidth', 0.5)
xlabel('time steps (k)'),ylabel('x_i(1)')
title('regulated state synchronization, case III')
legend('x_1', 'x_2', 'x_3', 'x_r')

subplot(4,1,2);
plot(k,x2, 'LineWidth', 1)
hold on
plot(k,x_r(2,:), 'LineWidth', 0.5)
xlabel('time steps (k)'),ylabel('x(2)'), ylabel('x_i(2)')
legend('x_1', 'x_2', 'x_3', 'x_r')

subplot(4,1,3);
plot(k,x3, 'LineWidth', 1)
hold on
plot(k,x_r(3,:), 'LineWidth', 0.5)
xlabel('time steps (k)'),ylabel('x_i(3)')
legend('x_1', 'x_2', 'x_3', 'x_r')

subplot(4,1,4);
plot(k,x - x_R, 'LineWidth', 1)
xlabel('time steps (k)'),ylabel('error')

