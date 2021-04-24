clc
clear all
close all


flag=0;
ID_im=5;
for i=1 : ID_im
    flag=1;
    switch i
        case 1
            img_orginal = imread('01.jpg') ;
            img_orginal = double(rgb2gray(img_orginal));
            c1 = -2;
            phi = ones(size(img_orginal)).*c1;
            phi(22:52,43:53) = -c1;
            figure()
            subplot(1,3,1)
            imshow(img_orginal,[])
            hold on
            contour(phi,[0 0],'r','LineWidth',2)
            itr   = 4500;
            landa = 0.99995;
            sigma = 11;
            phi   = proposed_method( img_orginal,phi,itr,landa,sigma );
            
        case 2
            img_orginal = imread('02.jpg') ;
            img_orginal = double(rgb2gray(img_orginal));
            c1  = -2;
            phi = ones(size(img_orginal)).*c1;
            phi(55:65,50:55) = -c1;
            itr = 1700;
            figure()
            subplot(1,3,1)
            imshow(img_orginal,[])
            hold on
            contour(phi,[0 0],'r','LineWidth',2)
            landa = 0.9995;
            sigma = 4;
            phi = proposed_method( img_orginal,phi,itr,landa,sigma );
        case 3
            img_orginal = imread('03.jpg') ;
            img_orginal = double(rgb2gray(img_orginal));
            c1  = -2;
            phi = ones(size(img_orginal)).*c1;
            phi(42:52,55:65) = -c1;
            itr = 12000;
            figure()
            subplot(1,3,1)
            imshow(img_orginal,[])
            hold on
            contour(phi,[0 0],'r','LineWidth',2)
            landa = 1;
            sigma = 4;
            phi = proposed_method( img_orginal,phi,itr,landa,sigma );

            
        case 4
            img_orginal = imread('04.jpg');
            img_orginal = double(rgb2gray(img_orginal));
            c1=-2;
            phi1 = ones(size(img_orginal)).*c1;
            phi1(5:70,5:115) = -c1;
            phi2 = ones(size(img_orginal)).*c1;
            phi2(20:60,35:75) = -c1;
            
            figure()
            subplot(1,3,1)
            imshow(img_orginal,[])
            hold on
            contour(phi1,[0 0],'b','LineWidth',2)
            contour(phi2,[0 0],'r','LineWidth',2)
            
            sigma = 14;
            landa = 0.9995;
            meu = 0.1;
            itr = 600;
            [phi1,phi2] = four_phase_formulation(img_orginal,phi1,phi2,itr,meu,landa,sigma);
            flag=1;
            
        case 5
            
           img_orginal = imread('03.jpg') ;
            itr=300;
            sigma=4;
            phi = selective_segmentation_model(img_orginal,itr,sigma);
            
    end
    
    if i==4
        subplot(1,3,3)
        imshow(img_orginal,[])
        hold on
        contour((phi1(:,:,1)),[0 0],'b','LineWidth',2)
        contour((phi2(:,:,1)),[0 0],'r','LineWidth',2)
    end
        
        if i==5
            subplot(1,2,2)
            imshow(img_orginal,[])
            hold on
            contour((phi(:,:,1)),[0 0],'b','LineWidth',2)
        end
            
        if i==1 | i==2 | i==3
            
            subplot(1,3,3)
            imshow(img_orginal,[])
            hold on
            contour((phi(:,:,1)),[0 0],'b','LineWidth',2)  
       
        end
end
