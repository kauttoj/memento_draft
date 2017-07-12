load hand
% display data
figure;
plot(p(:,1,1),p(:,2,1),'LineWidth',4)
axis('equal');




% compute the mean and PCs of these 9 hands
 

[mean_shape pc_shape std_shape pc_projection new_pt_out]=tangent_pca_shape(p);
 
figure
pareto(std_shape.^2)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

pc_index=2;

figure;
for i=-3:0.5:3
    temp = cos(i*std_shape(pc_index))*mean_shape+sin(i*std_shape(pc_index))*pc_shape(:,:,pc_index);
    plot(temp(:,1),temp(:,2),'LineWidth',4)
    pause(0.9)
end

    


