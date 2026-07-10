function X = proj_oblique(Y)
    [m,r]=size(Y);
    X=zeros(m,r);

    for i=1:r
        Yi=Y(:,i);
        X(:,i)=proj_unit_ball(Yi);
    end
end

function x = proj_unit_ball(y)
    if any(y,'all')
        x=y/norm(y);
    else
        n=length(y);
        r=rand(n,1);
        x=r/norm(r);
    end    
end
