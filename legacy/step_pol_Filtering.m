% I :  given image
% grammes: size, 
% DTheta: angle vector
% DPaxos:  vector of width of main wave
% typeOfFilter : type of step filter in {0,1} : 1 step function, 0 Polynomial filter 
% Y: filter output 



function [Y]= step_pol_Filtering(I,grammes, DTheta,DPaxos,typeOfFilter)  

gabors = length(DTheta);
gabors2 = length(DPaxos);



for j=1:gabors2,
    bw = DPaxos(j);
   
    for i=1:gabors,
        theta = DTheta(i);
        if typeOfFilter == 1,
            gb=mask_fn(grammes,bw,theta,0);
        else
            gb=mask_fn_smooth(grammes,bw,theta,0);
        end
        gb = gb / sum(sum(gb.^2)); %normalazation in energy

        G = conv2(I,gb,'same'); %idio size me I
        if i == 1 && j == 1,
            Y = abs(G);
        else
            Y = max(abs(Y),abs(G));
        end
    end
end

del = 8;
Y(1:del,:) = 0;
Y(:,1:del) = 0;
Y(size(I,1)-[0:del-1],:) = 0;
Y(:,size(I,2)-[0:del-1]) = 0;


i = 1;
if typeOfFilter == 1,
    gb=mask_fn(grammes,bw,pi/4,1);
else
    gb=mask_fn_smooth(grammes,bw,pi/4,1);
end


function gb=mask_fn(grammes,bw,theta,toPlot)
% bw    = bandwidth, (1)
% gamma = aspect ratio, (0.5)
%
% lambda= wave length, (>=2)
% theta = angle in rad, [0 pi)
 

sz=grammes;
if mod(sz,2)==0, sz=sz+1;end


[x y]=meshgrid(-fix(sz/2):1:fix(sz/2),fix(sz/2):-1:fix(-sz/2));


% Rotation 
x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);
 
for i=1:length(x),
    for j=1:length(y),
        d = abs(x_theta(i,j));
        u = x_theta(i,j);
        v = y_theta(i,j);
        gb(i,j) = 0;
        if d < bw,
            gb(i,j) = 1;
        elseif d <= 3*bw
            gb(i,j) = -1;
        end
    end
end

[n1 n2] = find(gb == -1);
[p1 p2] = find(gb == 1);


nv_n1 = -length(p1)/length(n1);

for i=1:length(n1),
    gb(n1(i),n2(i))= nv_n1;
end

gb = gb/sum(sum(gb.^2));

%gb=exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);

if toPlot == 1,
 
 figure;
 imagesc(gb);
 colormap(gray);
 title('Filter');

mean(mean(gb))
end


function gb=mask_fn_smooth(grammes,bw,theta,toPlot)
% bw    = bandwidth, (1)
% grammes = filter size
%
% toPlot= 1 for plotting  
% theta = angle in rad, [0 pi)
 

sz=grammes;
if mod(sz,2)==0, sz=sz+1;end


[x y]=meshgrid(-fix(sz/2):1:fix(sz/2),fix(sz/2):-1:fix(-sz/2));


% Rotation 
x_theta=x*cos(theta)+y*sin(theta);
y_theta=-x*sin(theta)+y*cos(theta);

    
for i=1:length(x),
    for j=1:length(y),
        d = abs(x_theta(i,j));
        u = x_theta(i,j);
        v = y_theta(i,j);
        gb(i,j) = 0;
        if d < bw,
            gb(i,j) = 1-((d/bw)^2);
        elseif d <= 3*bw
            gb(i,j) = 0.5*(-1+(((d/bw)-2)^2));
        end
    end
end


[u1 v1] = find(gb > 0);
[u2 v2] = find(gb < 0);

m1 = 0;
N1 = length(u1);
for i=1:length(u1),
    m1 = m1+gb(u1(i),v1(i));
end
m1 = m1/N1;

m2 = 0;
N2 = length(u2);
for i=1:length(u2),
    m2 = m2+gb(u2(i),v2(i));
end
m2 = m2/N2;

c = -m1*N1/(m2*N2);

for i=1:length(u2),
    gb(u2(i),v2(i)) = c*gb(u2(i),v2(i));
end

gb = gb/sum(sum(gb.^2));

%gb=exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi/lambda*x_theta+psi);

if toPlot == 1,
 
 figure;
 imagesc(gb);
 colormap(gray);
 title('Filter');

figure;
 mesh(gb);
  title('Filter in 3D view');
  
mean(mean(gb))
end


