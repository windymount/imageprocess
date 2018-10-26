function oimage = tv(input,kernel,mu,delta,rounds,lambda,ori)
%mu coefficient for augmented Lagrangian
input = im2double(input);
    n = size(input);
    u = zeros(n);
    z = zeros(n(1),2*n(2));
    L = [0 1 0;1 -4 1;0 1 0];
    G1=[0 0 0;0 1 0;0 -1 0];
    G2= [0 0 0;-1 1 0;0 0 0];
    y = [imfilter(input,G1,"conv","circular") imfilter(input,G2,"conv","circular")];
    for i = 1:rounds
        temp = imfilter(input,kernel,"conv","circular");
        g = y-z;
        temp = temp - mu*(imfilter(g(:,1:n(2)),G2',"conv","circular")+imfilter(g(:,n(2)+1:2*n(2)),G1',"conv","circular"));
        f =@(x)(imfilter(imfilter(x,kernel,"conv","circular"),kernel,"conv","circular")-mu*imfilter(x,L,"conv","circular"));
        u = cg(f,temp,400,0.001,n);
        disp(psnr(u,ori))
        temp = [imfilter(u,G1,"conv","circular") imfilter(u,G2,"conv","circular")]+z;
        v = norm(temp);
        y = temp/v*max(0,v-lambda/mu);
        z = z+delta*(temp-z-y);
    end
    oimage = u;
    
