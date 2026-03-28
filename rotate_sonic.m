function [A] = rotate_sonic(A)


for i = 1:size(A.U,2)
    u1 = A.U(:,i);
    v1 = A.V(:,i);
    w1 = A.W(:,i);
    [u2,v2] = rot_atan2(u1, v1);
    [u3,w2] = rot2_atan2(u2, w1);


    A.U(:,i) = u3; 
    A.V(:,i) = v2; 
    A.W(:,i) = w2; 

end

end



function [u,w]=rot2_atan2(utmp,wtmp)
w1mean = nanmean(wtmp);
u1mean = nanmean(utmp);

phi = atan2(w1mean, u1mean);

u = utmp * cos(phi) + wtmp * sin(phi);
w = -utmp * sin(phi) + wtmp * cos(phi);
end

function [u,v]=rot_atan2(u1,v1)
u=0*u1;v=0*v1;dum=0*u1;
v1mean=nanmean(v1);
u1mean=nanmean(u1);
theta=atan2(v1mean,u1mean);
%ws=sqrt((nanmean(v1))^2+(nanmean(u1))^2);
costheta=cos(theta);
sintheta=sin(theta);

u=u1*costheta+v1*sintheta;
v=-u1*sintheta+v1*costheta;
end

