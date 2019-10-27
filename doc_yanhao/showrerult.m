Tlc = ...
[  -0.0147587         -0   0.999891 -0.0481065;
 -0.999891         -0 -0.0147587  -0.016157;
         0         -1          0     0.0574;
         0          0          0          1;...
]


% x1 = [0.530 -0.191 0.000;
% -0.034 -3.346 0.000;...
% ]
% x2 = [
%     0.499 -0.364 0.000;
% 0.663 0.553 0.000;
% ]
% 
% 
% plot(x1(:,1),x1(:,2),'b-'); hold on
% plot(x2(:,1),x2(:,2),'r-') 


showresult(Tlc)

function [] = showresult(Tlc)


p_camer = Tlc(1:3,4);
R_camer = Tlc(1:3,1:3);

p_laser = [0;0;0];
R_laser = eye(3);


figure
axisl = 0.05;

% plot laser
xdirCamer = p_camer+axisl*R_camer(:, 1);
xdirCamer = [p_camer xdirCamer];
ydirCamer = p_camer+axisl*R_camer(:, 2);
ydirCamer = [p_camer ydirCamer];
zdirCamer = p_camer+axisl*R_camer(:, 3);
zdirCamer = [p_camer zdirCamer];

plot3(xdirCamer(1,:), xdirCamer(2,:), xdirCamer(3,:), 'Color', 'red'); hold on;
plot3(ydirCamer(1,:), ydirCamer(2,:), ydirCamer(3,:), 'Color', 'green'); hold on;
plot3(zdirCamer(1,:), zdirCamer(2,:), zdirCamer(3,:), 'Color', 'blue'); hold on;
 
% plot camera
xdirLaser = p_laser+axisl*R_laser(:, 1);
xdirLaser = [p_laser xdirLaser];
ydirLaser = p_laser+axisl*R_laser(:, 2);
ydirLaser = [p_laser ydirLaser];
zdirLaser = p_laser+axisl*R_laser(:, 3);
zdirLaser = [p_laser zdirLaser];

plot3(xdirLaser(1,:), xdirLaser(2,:), xdirLaser(3,:), 'Color', 'red'); hold on;
plot3(ydirLaser(1,:), ydirLaser(2,:), ydirLaser(3,:), 'Color', 'green'); hold on;
plot3(zdirLaser(1,:), zdirLaser(2,:), zdirLaser(3,:), 'Color', 'blue'); hold on; 

view(3)

text(p_camer(1)+0.005, p_camer(2)+0.005, p_camer(3)+0.005,'camer');
text(p_laser(1)+0.005,p_laser(1)+0.005,p_laser(1)+0.005,  'laser');


axis equal


end
