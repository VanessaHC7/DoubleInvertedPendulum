function Plot_Inverted_Pendulum(t,q,N,l,obj)

% Rotation matrix function  

%R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
R = @(theta) [cos(pi/2-theta) -sin(pi/2-theta); sin(pi/2-theta) cos(pi/2-theta)];


npoints = 5;

lv(1,:) = -1:2/(npoints-1):1;
for i = 2:(N+1)
    lv(i,:) = 0:l(i)/(npoints-1):l(i);
end

for j = 1:npoints
    x(:,1,j) = q(:,1) + lv(1,j);
    y(:,1,j) = zeros(size(q(:,1)));
end

for i = 2:(N+1)
    for k = 1:length(t)
        for j = 1:npoints
            x(k,i,j) = x(k,i-1,(npoints+1)/2) + [1 0]*R(q(k,i-1))*[l(i-1)/2 0]' + [1 0]*R(q(k,i))*[lv(i,j) 0]';
            y(k,i,j) = y(k,i-1,(npoints+1)/2) + [0 1]*R(q(k,i-1))*[l(i-1)/2 0]' + [0 1]*R(q(k,i))*[lv(i,j) 0]';
        end
    end
end


obj.Quality = 100;
obj.FrameRate = 24;
open(obj)

for i = 1:1e0:length(t)
    hold on
    %plot(x(1:i,1:end,(1+npoints)/2),y(1:i,1:end,(1+npoints)/2),'linewidth',1)
    for k = 1:(N+1)
        plot(squeeze(x(i,k,:)),squeeze(y(i,k,:)),'-','linewidth',4)
    end
    xlim([(min(x(i,1,:))-5) (max(x(i,1,:))+5)])
    ylim([-6 6])
    grid on
    pause(0.0001)

    f = getframe(gcf);
    writeVideo(obj,f);

    if i ~= length(t)
        clf
    end
end

close(obj)

end