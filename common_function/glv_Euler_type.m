function [dx,dxx]=glv_Euler_type(x,A,r,time,Type,h1,h2)
N = size(A,1);

nt = length(time);
dt = time(2) - time(1);

dxx = x;
dx = x;

switch Type
    case 1  % Functional type 1
        for i = 2:nt;
            dxx = dxx + (A'* dxx + r').*dxx * dt;
            dx(:,i) = dxx;
        end
    case 2 % Functional type 2
        for i = 2:nt;
            dxx = dxx + (A'* (dxx./(1+h1*dxx)) + r').*dxx * dt;
            dx(:,i) = dxx;
        end
    case 3 % Functional DeAngelis-Beddington (DB)
        for i = 2:nt;
            for j = 1 : N
                dxx(j,1) = dxx(j,1) + (A(:,j)'* (dxx./(1+h1*dxx+h2*dxx(j,1))) + r(j)).*dxx(j,1) * dt;
            end
            dx(:,i) = dxx;
        end
    case 4 % Functional Crowley-Martin (CM)
        for i = 2:nt;
            for j = 1 : N
                dxx(j,1) = dxx(j,1) + (A(:,j)'* (dxx./((1+h1*dxx)*(1+h2*dxx(j,1)))) + r(j)).*dxx(j,1) * dt;
            end
            dx(:,i) = dxx;
        end
end
dx = dx';
end