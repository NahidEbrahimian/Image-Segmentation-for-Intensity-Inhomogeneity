function  [phi1,phi2] = four_phase_formulation(Z0,phi1,phi2,iteration,meu,landa,sigma)

w = abs(4*sigma+1);
kernel = fspecial('gaussian', w,sigma);

[Row1, Column1] = size(phi1);
[Row2, Column2] = size(phi2);
epsilon = 1;

subplot(1,3,2)
imshow(Z0,[])


for itr =1:iteration
    
    hold on
    %     for i=1 : 4
    %         switch i
    %             case 1
    %                 [x1,y1] = find(phi1>0);
    %                 [x2,y2] = find(phi2>0);
    %                 R=findd(x1,y1,x2,y2);
    %                 R11=R;
    %             case 2
    %                 [x1,y1] = find(phi1>0);
    %                 [x2,y2] = find(phi2<0);
    %                 R=findd(x1,y1,x2,y2);
    %                 R01=R;
    %             case 3
    %                 [x1,y1] = find(phi1<0);
    %                 [x2,y2] = find(phi2>0);
    %                 R=findd(x1,y1,x2,y2);
    %                 R10=R;
    %             case 4
    %                 [x1,y1] = find(phi1<0);
    %                 [x2,y2] = find(phi2<0);
    %                 R=findd(x1,y1,x2,y2);
    %                 R11=R;
    %
    %         end
    %     end
    
    P = phi1;
    P([1 Row1],[1 Column1]) = P([3 Row1-2],[3 Column1-2]);
    P([1 Row1],2:end-1) = P([3 Row1-2],2:end-1);
    P(2:end-1,[1 Column1]) = P(2:end-1,[3 Column1-2]);
    phi1 = P;
    
    P = phi2;
    P([1 Row2],[1 Column2]) = P([3 Row2-2],[3 Column2-2]);
    P([1 Row2],2:end-1) = P([3 Row2-2],2:end-1);
    P(2:end-1,[1 Column2]) = P(2:end-1,[3 Column2-2]);
    phi2 = P;
    
    
    Delta_eps1 = (epsilon/pi)./(epsilon^2 + phi1.^2);
    epsilon = 1;
    H_eps1 = 0.5*(1 + (2/pi)*atan(phi1 ./ epsilon));
    
    Delta_eps2 = (epsilon/pi)./(epsilon^2 + phi2.^2);
    epsilon = 1;
    H_eps2 = 0.5*(1 + (2/pi)*atan(phi2 ./ epsilon));
    
    
    M1=H_eps1.*H_eps2;
    M2=H_eps1.*(1-H_eps2);
    M3=(1-H_eps1).*H_eps2;
    M4=(1-H_eps1).*(1-H_eps2);
    
    
    
    F1 = conv2(Z0 .* M1,kernel,'same')./ conv2(M1,kernel,'same');
    F2 = conv2(Z0 .* M2,kernel,'same')./ conv2(M2,kernel,'same');
    F3 = conv2(Z0 .* M3,kernel,'same')./ conv2(M3,kernel,'same');
    F4 = conv2(Z0 .* M4,kernel,'same')./ conv2(M4,kernel,'same');
    
    Z  = F1 .* M1 + F2 .* M2 + F3 .* M3 + F4 .* M4;
    
    D1  = ((F1-F2-F3+F4).*H_eps2)+F2-F4;
    D2  = ((F1-F2-F3+F4).*H_eps1)+F3-F4;
    
    [G_x_phi1,G_y_phi1] = gradient(phi1);
    xx = sqrt(G_x_phi1 .^ 2 + G_y_phi1 .^ 2+eps);
    G_x_phi1 = G_x_phi1./xx;
    G_y_phi1 = G_y_phi1./xx;
    kapa1 = divergence(G_x_phi1,G_y_phi1);
    
    [G_x_phi2,G_y_phi2] = gradient(phi2);
    xx = sqrt(G_x_phi2 .^ 2 + G_y_phi2 .^ 2+eps);
    G_x_phi2 = G_x_phi2./xx;
    G_y_phi2 = G_y_phi2./xx;
    kapa2 = divergence(G_x_phi2,G_y_phi2);
    
    temp1  = D2 .* (((2*(1-landa)) .* (Z - double(Z0))) + (landa.* (Z - double(Z0)) ./ (Z .^ 2)));
    phi2 = phi2 + Delta_eps2 .*(meu .* kapa2) - temp1;
    
    temp  = D1 .* (((2*(1-landa)) .* (Z - double(Z0))) + (landa.* (Z - double(Z0)) ./ (Z .^ 2)));
    phi1 = phi1 + Delta_eps1 .*(meu .* kapa1) - temp;
    
    if(rem(itr,50)==0)
     contour(phi1,[0 0],'b')
      contour(phi2,[0 0],'r')
    end
end




end

