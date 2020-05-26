function [dx,dxx]=glv_Euler_time(x,A,r,time)
N = size(A,1);

nt = length(time);
dt = time(2) - time(1);

dxx = x;
dx = x;
for i = 2:nt;
    dxx = dxx + (A'* dxx + r').*dxx * dt;
    dx(:,i) = dxx;
end
end