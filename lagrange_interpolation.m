function yy=lagrange_interpolation(X,pointX,pointY,N)
if (length(pointY)~=length(pointX))
    fprintf(1,'\nERROR!\nPOINTX and POINTY must have the same number of elements\n');
    y=NaN;
else
    X=X(:);
    F=1:length(pointY);
    F=F(:);
    for jj=1:length(X)
        t=[abs(pointX(:)-X(jj)),F];
        t=sortrows(t);
        t=t(1:N,2);
        pointx = pointX(t);
        pointy = pointY(t);
        x=X(jj);
        n=length(pointx);
        L=ones(n,length(x));
        for i=1:n
            for j=1:n
                if (i~=j)
                    L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
                end
            end
        end
        y=0;
        for i=1:n
            y=y+pointy(i)*L(i,:);
        end
        yy(jj,1)=y;
    end
end
end