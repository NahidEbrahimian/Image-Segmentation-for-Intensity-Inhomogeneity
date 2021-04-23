function  phi = selective_segmentation_model(Z0,iteration,sigma)

epsilon=1;

w = abs(4*sigma+1);
kernel = fspecial('gaussian', w,sigma);


pts = readPoints(Z0, 13);
x=round(pts(1,:));
y=round(pts(2,:));
x(1,11)=x(1);
y(1,11)=y(1);

Z0 = rgb2gray(Z0);
Z0 = double(Z0);

c1=-2;
polygon= zeros(size(Z0)).*c1;
polygon=roipoly(polygon,x,y);

V1=bwarea(polygon);
[l,w]=size(polygon);
S=l*w;
V2=S-V1;

c1=-2;
phi= ones(size(Z0)).*c1;
phi=roipoly(polygon,x,y);
% phi=phi.*4;
% phi=phi-2;
phi=bwdist(~phi)-bwdist(phi)-phi;


Delta_eps = (epsilon/pi) ./ (epsilon^2 + phi.^2);
epsilon = 1;
H_eps = 0.5 * (1 + (2/pi) * atan(phi ./ epsilon));

landa1=1;
landa2=1;
landa=0.9;
[len,wit]=size(Z0);
meu=0.0005*len*wit;
v=0.03;
% tao=1;
%
%
%
% dx = 0;
% for i = 1 : n1
%     a  = 1 - e^(-(((x-x(i))^2) / 2*(tao^2)))
%     dx = dx*a;
% end

dx=phi;

[x,y] = gradient(phi);
xx = sqrt(x .^ 2 + y .^ 2+eps);
g = 1/1+xx;
w = dx.*g;

x1 = (w.*x) ./ xx;
y1 = (w.*y) ./ xx;
divergen = divergence(x1,y1);


c1 = conv2(Z0 .* H_eps,kernel,'same')./ conv2(H_eps,kernel,'same');
c2 = conv2(Z0 .* (1 - H_eps),kernel,'same')./ conv2((1 - H_eps),kernel,'same');

F1 = conv2((Z0.*H_eps),kernel,'same')./ conv2(H_eps,kernel,'same');
F2 = conv2(Z0 .* (1 - H_eps),kernel,'same')./ conv2((1 - H_eps),kernel,'same');

D  = F1 - F2;
Z  = F1 .* H_eps + F2 .* (1-H_eps);

temp = meu .* divergen - (landa1 .* (Z0-c1).^2 - landa2 .* (Z0-c2).^2) -  v.* ((conv2((dx .* H_eps),kernel,'same'))-V1 - (conv2((dx .* (1-H_eps)),kernel,'same'))-V2);


phi = phi + Delta_eps .* ( temp - D .* (((2*(1-landa)) .* (Z - double(Z0))) + (landa.* (Z - double(Z0)) ./ (Z .^ 2)))) ;


end

