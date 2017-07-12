clc;
clear all;
close all;

N=80;
TR=1.56;

t1 = 20:TR:1000 + rand();
shift = 87 + 3*rand();
t2 = t1 + shift;

t1 = t1(1:N);
t2 = t2(1:N);

ind1 = round(N/5 - 3:(N-N/5 + 3)) - 10;
ind2 = ind1 + 20;

t1_new = t1(1)+(t1(end)-t1(1))/5 - 5;
t2_new = t2(1)+(t2(end)-t2(1))/5 + 5;

t1_new = t1_new:(TR/2):( t1_new + (t1(end)-t1(1))*3/5);
t2_new = t2_new:(TR/2):( t2_new + (t2(end)-t2(1))*3/5);

y1 = detrend(randn(1,N));
y2 = detrend(randn(1,N));
y2(ind2) = y1(ind1)*0.75 + 0.25*randn(1,length(ind2));

y1_new = interp1(t1,y1,t1_new,'linear');
y2_new = interp1(t2,y2,t2_new,'linear');

figure;
subplot(2,1,1);hold on
plot(t1,y1,'bo-');
plot(t1_new,y1_new,'rx-');
title('Nonlinear Memento')
legend('Original (rate=TR)','Interpolated (rate=2TR)')
xlabel('Time [s]')
axis tight;
box on

subplot(2,1,2);hold on
plot(t2,y2,'bo-');
plot(t2_new,y2_new,'rx-');
title('Linear Memento')
legend('Original (rate=TR)','Interpolated (rate=2TR)')
xlabel('Time [s]')
axis tight;
box on