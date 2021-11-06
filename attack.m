function f1=attack()
for i=16:16:64
    f1=imread(strcat('C:\Users\USER\Desktop\',int2str(1),'.jpg'));
    gau=imnoise(f1,'gaussian',0,.002);
   
    spk=imnoise(f1,'speckle',0.01);
    
    salt=imnoise(f1,'salt & pepper',0.5);
  
    po=imnoise(f1,'poisson');
    
    
    [row,col]=size(f1);
    for j=200:400
            f1(i,j)=255;
        end
end
    q=90;
    tmp=int2str(q);
    name=strcat('cmprs_',tmp,'%.jpeg');
    imwrite(f1,name,'jpeg','Bitdepth',8,'quality',q);
    f2=imread(name);
    S=decorrstretch(f1,'tol',.02);
    p=wiener2(f1,[2,2]);
    watermarked_image_int=uint8(p);
    imwrite(watermarked_image_int,strcat('C:\Users\USER\Desktop\',int2str(i),'-2.jpg'),'jpg');
    disp(i);
end