
function ext= extraction_final()
for k=16:16:64
    blocksize=8;
    pn_sequence_search='T';
    lowband=[ 1,1 1,0;
              1,1,0,0;
              1,0,0,0;
              0,0,0,0];
    h=[0 1 0; 1 -4 1; 0 1 0];
    file_name=strcat(strcat('C:\Users\USER\Desktop\64-2.jpg'));
    watermarked_image=double(imread(file_name));
    I=imfilter(watermarked_image,h);
    h1=fspecial('log');
    watermarked_image=imfilter(I,h1);
    
    Mw=size(watermarked_image,1);
    Nw=size(watermarked_image,2);
   
imagesc(watermarked_image, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
disp( ' ' );
nlevs = [0, 2] ;        % Decomposition level
pfilt = 'pkva' ;              % Pyramidal filter
dfilt = 'pkva' ;
A=size(watermarked_image);
disp(A);
if length(nlevs) == 0
    y = {watermarked_image};
    
else
    % Get the pyramidal filters from the filter name
    [h, g] = pfilters(pfilt);
    
    if nlevs(end) ~= 0
        % Laplacian decomposition
        [xlo, xhi] = lpdec(double(watermarked_image), h, g);
    
        
        switch dfilt        % Decide the method based on the filter name
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}   
                % Use the ladder structure (whihc is much more efficient)
                xhi_dir = dfbdec_l(xhi, dfilt, nlevs(end));
            
            otherwise       
                % General case
                xhi_dir = dfbdec(xhi, dfilt, nlevs(end));                
        end
        
    else        
        
        [xlo, xLH, xHL, xHH] = wfb2dec(watermarked_image, h, g);
        xhi_dir = {xLH, xHL, xHH};
    end
    
    % Recursive call on the low band
    ylo = pdfbdec(double(xlo), pfilt, dfilt, nlevs(1:end-1));
end
    data=cell2mat([ylo{:}]);
cov_obj3=imresize(data, [512, 512]);
max_message=Mw*Nw/(blocksize^2);

file_name='C:\Users\USER\Desktop\peppers.png';
orig_watermark1=round((double(imread(file_name)))./256);

Mo=size(orig_watermark1,1);
No=size(orig_watermark1,2);
rand('state',16);
message_vector=ones(1,Mo*No);
    pn_sequence_one=round(2*(rand(1,sum(sum(lowband)))-0.5));
    pn_sequence_zero=round(2*(rand(1,sum(sum(lowband)))-0.5));
    if(pn_sequence_search=='T')
        while(corr2(pn_sequence_one,pn_sequence_zero)>-0.55)
         pn_sequence_one=round(2*(rand(1,sum(sum(lowband)))-0.5));
         pn_sequence_zero=round(2*(rand(1,sum(sum(lowband)))-0.5));
        end
    end
    A=imresize(pn_sequence_one,[1 40]);
    disp(A);
     B=imresize(pn_sequence_zero,[1 40]);
    disp(B);
    n=256;
    M = rescale(cov_obj3,n);
options.null = 0;
options.finest = 1;
options.nbscales = 4;
options.nbangles_coarse = 16;
options.is_real = 1;
options.n = n;
    MW = perform_curvelet_transform(M, options);
    cov_obj4=cell2mat([MW{1}]);
    
    for(kk=1:length(message_vector))
        
ii=1;
jj=1;
ll=1;
% k is the value of alpha
blocksize=8;

for ii=1:blocksize
    for jj=1:blocksize
        if lowband(jj)==1
               sequence(ll)=cov_obj4(jj,ii);
                ll=ll+1;
    end
    end
    
    end
    end
S=size(sequence);
disp(S);
   corr_one(kk)=corr2(A,sequence);
   corr_zero(kk)=corr2(B,sequence);
    
if (corr_zero(kk)>=corr_one(kk))
    message_vector(kk)=0;
else
    message_vector(kk)=1;
end

             
messa=reshape(message_vector,Mo,No);
imagesc(messa, [0, 255]);
    [r,c] = size(messa);  rc = r*c;
    K = CHC_chaos(rc);
    C0 = 74;
    Pd = double(messa);  
  
    message= uint8(CHC_decrypt(K,C0,messa));
    cP_d = corr_diagonal(message);
    cP_h = corr_horizontal(message);
    cP_v = corr_vertical(message);
    imwrite(message,strcat('C:\Users\USER\Desktop\',int2str(k),'-w.jpg'),'jpg');
    figure(1);
    imshow(message,[]);
    title('Recovered Message');
    correlation=corr2(message,orig_watermark1)+1;
    disp(correlation);
end
end 
