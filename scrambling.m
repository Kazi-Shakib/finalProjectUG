clear
data{1,6} = [];
data{10,6} = [];
    im = imread('C:\Users\USER\Desktop\peppers.png'); 
    if ~exist('im', 'var')
    % barbara image: good for illustrating multiscale and directional
    % decomposition
    im = imread ('C:\Users\USER\Desktop\peppers.png') ;
end

% Show the input image
disp( 'Displaying the input image...');
clf;
imagesc(im, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
input( 'Press Enter key to continue...' ) ;
disp( ' ' );
    [r,c] = size(im);  rc = r*c;
    K = CHC_chaos(rc);
    C0 = 74;
    Pd = double(im);  
    Q = uint8(CHC_encrypt(Pd,K,C0));
    cQ_d = corr_diagonal(Q);
    cQ_h = corr_horizontal(Q);
    cQ_v = corr_vertical(Q);
    
    imshow(Q);
    P = uint8(CHC_decrypt(K,C0,double(Q)));
    cP_d = corr_diagonal(P);
    cP_h = corr_horizontal(P);
    cP_v = corr_vertical(P);
    imshow(P);