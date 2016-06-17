%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%simulator to simulate the trajectory for two robots and landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
clc
%%
%%%%%
%%simualtion space 
lx=1000;
ly=1000;
lz=2000;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generate the robot trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%generate the trajectory of robot 1 
%%%
%Starting point for R1
R1_start=[10;250;10];
%Ending point for R1
R1_end=[900;400;1700];
%distance traveled:
R1_distance=R1_end-R1_start;
%randomly devided the steps
%steps:
step_n=50;
R1_steps=zeros(3,step_n);
for i=1:3
    R1_steps(i,:)=rand(1,step_n);
    R1_steps(i,:)=R1_distance(i)/sum(R1_steps(i,:))*R1_steps(i,:);
end

% the trajectory for R1
R1_traj=zeros(3,step_n+1);
R1_traj(:,1)=R1_start;
for i=2:step_n+1
    R1_traj(:,i)=R1_traj(:,i-1)+R1_steps(:,i-1);
end

%%%%%%%%%%
%generate the trajectory for robot 2
%%%
%Starting point for R2
R2_start=[200;550;10];
%Ending point for R2
R2_end=[800;600;1900];
%distance traveled:
R2_distance=R2_end-R2_start;
%randomly devided the steps
%steps:
step_n2=50;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Generate the landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lm_N=600;
lm=zeros(3,lm_N);
lm(1,:)=lx*rand(1,lm_N);
lm(2,:)=ly*rand(1,lm_N);
lm(3,:)=lz*rand(1,lm_N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Generate the measurement between trajectory and landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sensor range for robot
R1_range=200;
R2_range=200;

R1_z=struct('ID',0);
for i=1:step_n+1
    R1_z(i).ID=0;
    for j=1:lm_N
        d=norm(R1_traj(:,i)-lm(:,j));
        if d<=R1_range && d>=20
            R1_z(i).ID=[R1_z(i).ID, j];
        end
    end
end

R2_z=struct('ID',0);
for i=1:step_n2+1
    R2_z(i).ID=0;
    for j=1:lm_N
        d=norm(R2_traj(:,i)-lm(:,j));
        if d<=R2_range && d>=20
            R2_z(i).ID=[R2_z(i).ID, j];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%draw the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot3(R1_traj(1,:),R1_traj(2,:),R1_traj(3,:),'r-');
hold on
plot3(R2_traj(1,:),R2_traj(2,:),R2_traj(3,:));
hold on
%plot3(lm(1,:),lm(2,:),lm(3,:),'k+');
for i=1:lm_N
    %text(lm(1,i)+20,lm(2,i)+20,lm(3,i)+20,num2str(i));
end
hold on
%%%draw connections between robot trajectory and landmark
%{
for i=1:step_n+1
    connection_num=size(R1_z(i).ID,2);
    if connection_num>1
        ID=R1_z(i).ID;
        for j=2:connection_num
            plot3([R1_traj(1,i) lm(1,ID(j))],[R1_traj(2,i) lm(2,ID(j))],[R1_traj(3,i) lm(3,ID(j))],'r-.');
            %text(lm(1,ID(j))+20,lm(2,ID(j))+20,lm(3,ID(j))+20,num2str(ID(j)));
            hold on
            
        end
    end
end

for i=1:step_n2+1
    connection_num=size(R2_z(i).ID,2);
    if connection_num>1
        ID=R2_z(i).ID;
        for j=2:connection_num
            plot3([R2_traj(1,i) lm(1,ID(j))],[R2_traj(2,i) lm(2,ID(j))],[R2_traj(3,i) lm(3,ID(j))],'--');
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





