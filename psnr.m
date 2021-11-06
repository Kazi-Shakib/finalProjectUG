function f=psnr()
for i=16:16:64
    f1=double(imread(strcat('C:\Users\USER\Desktop\64-2.jpg')));
    f2=double(imread(strcat('C:\Users\USER\Desktop\1.jpg')));
    [row,col]=size(f1);
    sum=0;
    for k=1:row
        for j=1:col
            sum=sum+(f2(k,j)-f1(k,j))^2;
        end
    end
    mse=(sum/(row*col));
    psnr=10* log10(double(255^2/mse));
    disp(psnr);
end 
end
    