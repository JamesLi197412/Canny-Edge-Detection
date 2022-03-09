clear asll;clc; close all;
% Task 1 : Implement Gaussian Blur
% Written by : Zhiyue

% Read an Image
filename ='test05.jpg';
Img = imread(filename);

% Original image
figure,imshow(Img);
title('original image');
I = double(Img);
rows = size(I,1);
columns = size(I,2);

% 5 X 5 Kernel for test
Kernel = 1/159 .* [2,4,5,4,2;4,9,12,9,4;5,12,15,12,5;4,9,12,9,4;2,4,5,4,2];

%Initialize
Output=zeros(size(I));
%Pad the vector with zeros
I = padarray(I,[2,2],0,'pre');

[OutputMatrix] = convolution(I,Kernel);
% Image after convolution
OutputMatrix_show = uint8(OutputMatrix);
figure,imshow(OutputMatrix_show);
title('Image After Convolution');


% Task 2
% Sobel kernel 
sobel_Gx = [-1,0,1;-2,0,2;-1,0,1];
sobel_Gy = [-1,-2,-1;0,0,0;1,2,1];

% Gx, Gy
Tempout = padarray(OutputMatrix,[1,1],0,'both');
[Gx] = convolution(Tempout,sobel_Gx);
[Gy] = convolution(Tempout,sobel_Gy);
figure,imshow(rescale(Gx, 0, 255), [0, 255]);
figure,imshow(rescale(Gy, 0, 255), [0, 255]);

% Task 3
G_mag = hypot(Gx,Gy);
% figure 3
figure,imshow(rescale(G_mag, 0, 255), [0, 255]);
title('Gradient magnitude');


rows = size(G_mag,1);
columns = size(G_mag,2);

% Compute the Edge Direction in degrees
G_orientation = zeros(rows,columns);
for i = 1:rows
    for j = 1:columns
        G_orientation(i,j) = atand(Gy(i,j)./Gx(i,j));
        if G_orientation(i,j)<0
            G_orientation(i,j) = abs(G_orientation(i,j));
        end
    end
end
        
rows_new = size(G_orientation,1);
columns_new = size(G_orientation,2);
% Round - up degree to 0,45,90,135
canny_edge = zeros(rows_new,columns_new);
for i=1:rows_new
    for j=1:columns_new
        if G_orientation(i,j) >337.5 || G_orientation(i,j) >=0 && G_orientation(i,j) < 22.5 || G_orientation(i,j) > 157.5 && G_orientation(i,j) <= 202.5
            canny_edge(i,j) = 0;
        end
        
        if (G_orientation(i,j) >= 22.5 && G_orientation(i,j) < 67.5) || (G_orientation(i,j) > 202.5 && G_orientation(i,j) <= 247.5)
            canny_edge(i,j) = 45;
        end
        
        if (G_orientation(i,j) >= 67.5 && G_orientation(i,j) < 112.5) || (G_orientation(i,j) > 247.5 && G_orientation(i,j) <= 292.5) 
            canny_edge(i,j) = 90;
        end 
        
        if (G_orientation(i,j) >= 112.5 && G_orientation(i,j) <= 157.5) || (G_orientation(i,j) > 292.5 && G_orientation(i,j) <= 337.5)
            canny_edge(i,j) = 135;
        end
    end
end

% figure 4
figure,imshow(canny_edge);
title('Canny edge detection');

rows = size(canny_edge,1);
columns = size(canny_edge,2);

% Non - maxima Suppression and thresholding
G_mark = G_mag;
G_mag = padarray(G_mag,[1,1],0,'both');
for i=2:rows+1
    for j=2:columns+1
        if canny_edge(i-1,j-1) == 0
            if (G_mag(i,j) < G_mag(i,j+1) || (G_mag(i,j) < G_mag(i,j-1)))
                G_mark(i-1,j-1) = 0;
            end
        end
        
        if canny_edge(i-1,j-1) == 45
            if (G_mag(i,j) < G_mag(i-1,j+1)) || (G_mag(i,j) < G_mag(i+1,j-1))
                G_mark(i-1,j-1) = 0;
            end
        end
        
        if canny_edge(i-1,j-1) == 90
            if ((G_mag(i,j) < G_mag(i-1,j)) || (G_mag(i,j) < G_mag(i+1,j)))
                G_mark(i-1,j-1) = 0;
            end
        end
        
         if canny_edge(i-1,j-1) == 135
            if (G_mag(i,j) < G_mag(i-1,j-1)) || (G_mag(i,j) < G_mag(i+1,j+1))
                G_mark(i-1,j-1) = 0;
            end
         end
    end
end

% Result after Non Maximum Suppression
figure,imshow(G_mark);
title('non - maxima Suppression');

% Threshold setting
% For Local Thresholding
thresh_high = 0.16*255;           
thresh_low = 0.13*255;              

X_Hyst = G_mark;

% Hysterisis
X_Hyst_Pad = padarray(X_Hyst,[1 1],0,'both');

for i=2:rows+1
    for j=2:columns+1
        if X_Hyst_Pad(i,j) >= thresh_high
            X_Hyst(i-1,j-1) = 1;
        end
        if X_Hyst_Pad(i,j) < thresh_high && X_Hyst_Pad(i,j) >= thresh_low
            if X_Hyst_Pad(i,j+1) >= thresh_high || X_Hyst_Pad(i-1,j) >=thresh_high || X_Hyst_Pad(i,j-1) >= thresh_high || X_Hyst_Pad(i+1,j) >= thresh_high || X_Hyst_Pad(i+1,j+1) >=thresh_high || X_Hyst_Pad(i-1,j+1) >= thresh_high || X_Hyst_Pad(i-1,j-1) >= thresh_high || X_Hyst_Pad(i+1,j-1) >= thresh_high
                X_Hyst(i-1,j-1) = 1;
            else
                X_Hyst(i-1,j-1) = 0;
            end  
        end
        if X_Hyst_Pad(i,j) < thresh_low
            X_Hyst(i-1,j-1) = 0;
        end
    end
end

image_final = X_Hyst;
figure,imshow(image_final);
title('Threshold');

