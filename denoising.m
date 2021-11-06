function f=psnr()
pfilt = '9-7';
dfilt = 'pkva';
nlevs = [0, 0, 4, 4, 5];    % Number of levels for DFB at each pyramidal level
th = 3;                     % lead to 3*sigma threshold denoising
rho = 3; 


im=imread(strcat('C:\Users\USER\Desktop\64-2.jpg'));
im = double(im) / 256;
y = pdfbdec(double(im), pfilt, dfilt, nlevs);
[c, s] = pdfb2vec(y);

% Threshold
% Require to estimate the noise standard deviation in the PDFB domain first 
% since PDFB is not an orthogonal transform
nvar = pdfb_nest(size(im,1), size(im, 2), pfilt, dfilt, nlevs);
sig = std(im(:));
sigma = sig / rho;
nim = im + sigma * randn(size(im));

cth = th * sigma * sqrt(nvar);
% cth = (4/3) * th * sigma * sqrt(nvar);
wth = th * sigma;
c = c .* (abs(c) > wth);

% Reconstruction
y = vec2pdfb(c, s);
wim = pdfbrec(y, pfilt, dfilt);

% Slightly different thresholds for the finest scale
fs = s(end, 1);
fssize = sum(prod(s(find(s(:, 1) == fs), 3:4), 2));
cth(end-fssize+1:end) = (4/3) * cth(end-fssize+1:end);

c = c .* (abs(c) > cth);

% Reconstruction
y = vec2pdfb(c, s);
cim = pdfbrec(y, pfilt, dfilt);

    [row,col]=size(im);
    sum1=0;
    for k=1:row
        for j=1:col
            sum1=sum1+(cim(k,j)-im(k,j))^2;
        end
    end
    mse=(sum1/(row*col));
    psnr=10* log10(double(255^2/mse));
    disp(psnr);

sum2=0;
    for k=1:row
        for j=1:col
            sum2=sum2+(wim(k,j)-im(k,j))^2;
        end
    end
    mse=(sum2/(row*col));
    psnr=10* log10(double(255^2/mse));
    disp(psnr);
 
imshow(cim);
    
