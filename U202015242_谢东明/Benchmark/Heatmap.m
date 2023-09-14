z = zeros(40,6);
for i = 1:40
    for j = 1:6
        z(i,j) = sum(Pred(:,i,j))/size(Pred(:,i,j),1);
    end
end
x = 1:40;
y = 1:6;
xx = 1:0.1:40;
yy = 1:0.05:6;
[X,Y] = meshgrid(x,y);
[XX,YY] = meshgrid(xx,yy);
Z = griddata(X,Y,z',XX,YY,'cubic');
for i = 1:size(Z,1)*size(Z,2)
    if Z(i) > 1
        Z(i) = 1;
    end
end
pcolor(XX,YY,Z)
xlabel('trial')
ylabel('block')
colorbar
shading flat