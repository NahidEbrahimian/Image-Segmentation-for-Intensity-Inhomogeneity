function  phi=proposed_method( Z0 ,phi,iteration,landa,sigma )

w = abs(4*sigma+1);
kernel = fspecial('gaussian', w,sigma);

[Row, Column] = size(phi);
epsilon = 1;
[len,wit]=size(Z0);
meu=0.000005*len*wit;
subplot(1,3,2)
imshow(Z0,[])

% subplot(1,3,3)
% imshow(Z0,[])
% hold on

for itr =1:iteration
    hold on
    P = phi;
    P([1 Row],[1 Column]) = P([3 Row-2],[3 Column-2]);
    P([1 Row],2:end-1) = P([3 Row-2],2:end-1);
    P(2:end-1,[1 Column]) = P(2:end-1,[3 Column-2]);
    phi = P;
    
    Delta_eps = (epsilon/pi)./(epsilon^2 + phi.^2);
    epsilon = 1;
    H_eps = 0.5*(1 + (2/pi)*atan(phi ./ epsilon));
    
    F1 = conv2((Z0.*H_eps),kernel,'same')./ conv2(H_eps,kernel,'same');
    F2 = conv2(Z0 .* (1 - H_eps),kernel,'same')./ conv2((1 - H_eps),kernel,'same');
    
    [x,y] = gradient(phi);
    xx = sqrt(x .^ 2 + y .^ 2+eps);
    x1 = x./xx;
    y1 = y./xx;
    kapa = divergence(x1,y1);
    
    D  = F1 - F2;
    Z  = F1 .* H_eps + F2 .* (1-H_eps);
    x  = D .* (((2*(1-landa)) .* (Z - double(Z0))) + (landa.* (Z - double(Z0)) ./ (Z .^ 2)));
    phi = phi + Delta_eps .*(meu .* kapa) - 1 / 3 .*(x);
    
    if rem(itr,300)==0

%         hold on        
        contour(phi,[0 0],'r')
    end
    
end

