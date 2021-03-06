%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%2 simulator to simulate the trajectory for two robots and landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
clc
%%
%%%%%
%%simualtion space 
lx=500;
ly=500;
lz=1200;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generate the robot trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%generate the trajectory of robot 1 and 2
%%%
%Starting point for R1
R1_start=[10;100;10];
R1_mid=[150;180;600];
%Ending point for R1
R1_end=[450;200;600];

%Starting point for R2
R2_start=[50;200;10];
%Ending point for R2
R2_end=[480;250;1000];



%distance traveled:
R1_dist1=R1_mid-R1_start;
R1_dist2=R1_end-R1_mid;
%randomly devided the steps
%steps:
step_n1=20;
step_n2=80;
step_n=step_n1+step_n2;
R1_steps=zeros(3,step_n);
for i=1:3
    R1_steps(i,1:step_n1)=rand(1,step_n1);
    R1_steps(i,step_n1+1:step_n)=rand(1,step_n2);
    R1_steps(i,1:step_n1)=R1_dist1(i)/sum(R1_steps(i,1:step_n1))*R1_steps(i,1:step_n1);
    R1_steps(i,step_n1+1:step_n)=R1_dist2(i)/sum(R1_steps(i,step_n1+1:step_n))*R1_steps(i,step_n1+1:step_n);
end



% the trajectory for R1
R1_traj=zeros(3,step_n+1);
R1_traj(:,1)=R1_start;
for i=2:step_n+1
    R1_traj(:,i)=R1_traj(:,i-1)+R1_steps(:,i-1);
end

R1_traj(3,:)=R1_traj(3,:)+5*randn(1,step_n+1);

%distance traveled:
R2_distance=R2_end-R2_start;
%randomly devided the steps
R2_steps=zeros(3,step_n);
for i=1:3
    R2_steps(i,:)=rand(1,step_n);
    R2_steps(i,:)=R2_distance(i)/sum(R2_steps(i,:))*R2_steps(i,:);
end

% the trajectory for R2
R2_traj=zeros(3,step_n+1);
R2_traj(:,1)=R2_start;
for i=2:step_n+1
    R2_traj(:,i)=R2_traj(:,i-1)+R2_steps(:,i-1);
end
%R2_traj(3,:)=R1_traj(3,:);
R2_traj(3,:)=R1_traj(3,:)+5*randn(1,step_n+1);

%%%%%%%%%%%%%%%%%%%%
%%fit the 3D curve
%%%%%%%%%%%%%%%%%%%%
x1=R1_traj(1,:)';
y1=R1_traj(2,:)';
z1=R1_traj(3,:)';
%T1=[x1 y1 x1.*y1 x1.*x1 y1.*y1 ones(step_n+1,1)]\z1;
%syms x1 y1
%z1=vpa([x1 y1 x1*y1 x1*x1 y1*y1 1]* T1, 4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generate the landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lm_N=100;
lm=zeros(3,lm_N);
lm(1,:)=lx*rand(1,lm_N);
lm(2,:)=ly*rand(1,lm_N);
lm(3,:)=100*rand(1,lm_N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generate the measurement between trajectory and landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sensor range for robot
R1_range=2000;
R2_range=2000;
vi_ang=pi/20;

step=20;


R1_z=struct('ID',0);
for i=step_n+1:-step:1
    R1_z(i).ID=0;
    for j=1:lm_N
        dr=norm(R1_traj(1:2,i)-lm(1:2,j));
        d=norm(R1_traj(:,i)-lm(:,j));
        dz=R1_traj(3,i)-lm(3,j);
        if dz>10 && dr<=dz*tan(vi_ang) && d<=R1_range
            R1_z(i).ID=[R1_z(i).ID, j];
        end
    end
end

R2_z=struct('ID',0);
for i=step_n+1:-step:1
    R2_z(i).ID=0;
    for j=1:lm_N
        dr=norm(R2_traj(1:2,i)-lm(1:2,j));
        d=norm(R2_traj(:,i)-lm(:,j));
        dz=R2_traj(3,i)-lm(3,j);
        if dz>10 && dr<=dz*tan(vi_ang) && d<=R2_range 
            R2_z(i).ID=[R2_z(i).ID, j];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%draw the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot3(R1_traj(1,:),R1_traj(2,:),R1_traj(3,:),'r-','LineWidth',3);
hold on
plot3(R2_traj(1,:),R2_traj(2,:),R2_traj(3,:),'LineWidth',3);
hold on
%plot3(lm(1,:),lm(2,:),lm(3,:),'k+');
for i=1:lm_N
    %text(lm(1,i)+20,lm(2,i)+20,lm(3,i)+20,num2str(i));
end
hold on
%%%draw connections between robot trajectory and landmark

for i=step_n+1:-step:1
    connection_num=size(R1_z(i).ID,2);
    if connection_num>1
        ID=R1_z(i).ID;
        for j=2:connection_num
            plot3([R1_traj(1,i) lm(1,ID(j))],[R1_traj(2,i) lm(2,ID(j))],[R1_traj(3,i) lm(3,ID(j))],'r-.*');
            %text(lm(1,ID(j))+20,lm(2,ID(j))+20,lm(3,ID(j))+20,num2str(ID(j)));
            hold on
            
        end
    end
end

for i=step_n+1:-step:1
    connection_num=size(R2_z(i).ID,2);
    if connection_num>1
        ID=R2_z(i).ID;
        for j=2:connection_num
            plot3([R2_traj(1,i) lm(1,ID(j))],[R2_traj(2,i) lm(2,ID(j))],[R2_traj(3,i) lm(3,ID(j))],'--o');
            %text(lm(1,ID(j))+20,lm(2,ID(j))+20,lm(3,ID(j))+20,num2str(ID(j)));
            hold on
        end
    end
end
%}

axis([0 lx 0 ly 0 lz]) 
grid on
xlabel('X');
ylabel('Y');
zlabel('Z');





