%%

function coeffs = cvt( cover_object, ~ )
if ~exist('cover_object', 'var')
    cover_object = imread ('C:\Users\USER\Desktop\barbara.png') ;
   
end


              % Directional filter


disp( 'Displaying the input image...');
clf;
imagesc(cover_object, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
disp( ' ' );
nlevs = [0, 2] ;        % Decomposition level
pfilt = 'pkva' ;              % Pyramidal filter
dfilt = 'pkva' ;
A=size(cover_object);
disp(A);
if length(nlevs) == 0
    y = {cover_object};
    
else
    % Get the pyramidal filters from the filter name
    [h, g] = pfilters(pfilt);
    
    if nlevs(end) ~= 0
        % Laplacian decomposition
        [xlo, xhi] = lpdec(double(cover_object), h, g);
    
        
        switch dfilt        % Decide the method based on the filter name
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}   
                % Use the ladder structure (whihc is much more efficient)
                xhi_dir = dfbdec_l(xhi, dfilt, nlevs(end));
            
            otherwise       
                % General case
                xhi_dir = dfbdec(xhi, dfilt, nlevs(end));                
        end
        
    else        
        
        [xlo, xLH, xHL, xHH] = wfb2dec(cover_object, h, g);
        xhi_dir = {xLH, xHL, xHH};
    end
    
    % Recursive call on the low band
    ylo = pdfbdec(double(xlo), pfilt, dfilt, nlevs(1:end-1));

    
end
% Display the coefficients
disp('Displaying the contourlet coefficients...') ;
imcoeff = showpdfb( ylo) ;
title('Contourlet coefficients');
input('Press Enter key to continue...' ) ;
disp(' ');

data=cell2mat([ylo{:}]);
   

cov_obj2=im2col(data,[8 8],'sliding');
cov_obj3=imresize(cov_obj2, [256, 256]);
A=size(cov_obj2);
disp(A);
imshow(cov_obj3);
% Contourlet decomposition
% Display the coefficients

n = 256;
name = 'barbara';
cover_object = imread ('C:\Users\USER\Desktop\barbara.png') ;
if ~exist('cover_object', 'var')
    % barbara image: good for illustrating multiscale and directional
    % decomposition
    cover_object = imread ('C:\Users\USER\Desktop\barbara.jpg') ;
end

% Show the input image
disp( 'Displaying the input image...');
clf;
imagesc(cover_object, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
input( 'Press Enter key to continue...' ) ;
disp( ' ' );
M = rescale(load_image(name,n));
options.null = 0;
options.finest = 1;
options.nbscales = 4;
options.nbangles_coarse = 16;
options.is_real = 1;
options.n = n;
MW = perform_curvelet_transform(M, options);
clf;

cov_obj4=cell2mat([MW{1}]);
cov_obj5=imresize(cov_obj4, [256, 256]);


%%
con=convn(cov_obj3,cov_obj5,'same');
imshow(con);
v=size(con);
disp(v);
%%
file_name='C:\Users\USER\Desktop\peppers.png';
mess=double(imread(file_name));

x=size(mess);
disp(x);
pn_sequence_search='T';
lowband=[ 1,1 ,1,0;
          1,1,0,0;
          1,0,0,0;
          0,0,0,0];
imagesc(mess, [0, 255]);
    [r,c] = size(mess);  rc = r*c;
    K = CHC_chaos(rc);
    C0 = 74;
    Pd = double(mess);  
    Q = uint8(CHC_encrypt(Pd,K,C0));
    cQ_d = corr_diagonal(Q);
    cQ_h = corr_horizontal(Q);
    cQ_v = corr_vertical(Q);
    

Mo=size(Q,1);
No=size(Q,2);


message=round(reshape(Q,Mo*No,1)./256);
messa = imresize(message, [64 64]);

%%
Mc=size(cover_object,1);
Nc=size(cover_object,2);
blocksize=8;
max_message=Mc*Nc/(blocksize^2);
message_vector=ones(1,max_message);
message_vector(1:length(message))=message;  

watermarked_image=cover_object;
rand('state',16);
    pn_sequence_one=round(2*(rand(1,sum(sum(lowband)))-0.5));
    pn_sequence_zero=round(2*(rand(1,sum(sum(lowband)))-0.5));
    if(pn_sequence_search=='T')
        while(corr2(pn_sequence_one,pn_sequence_zero)>-0.55)
         pn_sequence_one=round(2*(rand(1,sum(sum(lowband)))-0.5));
         pn_sequence_zero=round(2*(rand(1,sum(sum(lowband)))-0.5));
        end
    end
    %%Embedding Part

for(kk=1:length(message_vector))
ii=1;
jj=1;
ll=1;
k=0.5;
% k is the value of alpha
blocksize=8;
if(message_vector(kk)==0)
  
                con=con+k*pn_sequence_zero(ll);
                ll=ll+1;
        
   end
    
if (message_vector(kk)==1)
     
                con=con+k*pn_sequence_one(ll);
                ll=ll+1;
end
end    
T=0.01;
MWT=perform_thresholding(con,T,'hard');

Watermarked_image=perform_curvelet_transform(MWT,options);
clf;
s=size(Watermarked_image);
disp(s);

%inverse contourlet
function x = pdfbrec(Watermarked_image, pfilt, dfilt)
% PDFBREC   Pyramid Directional Filterbank Reconstruction
%
%	x = pdfbrec(y, pfilt, dfilt)
%
% Input:
%   y:	    a cell vector of length n+1, one for each layer of 
%       	subband images from DFB, y{1} is the low band image
%   pfilt:  filter name for the pyramid
%   dfilt:  filter name for the directional filter bank
%
% Output:
%   x:      reconstructed image
%
% See also: PFILTERS, DFILTERS, PDFBDEC

n = length(y) - 1;
if n <= 0
    x = Watermarked_image{1};
    
else
    % Recursive call to reconstruct the low band
    xlo = pdfbrec(Watermarked_image(1:end-1), pfilt, dfilt);
    
    % Get the pyramidal filters from the filter name
    [h, g] = pfilters(pfilt);
    
    % Process the detail subbands
    if length(Watermarked_image{end}) ~= 3
        % Reconstruct the bandpass image from DFB
        
        % Decide the method based on the filter name
        switch dfilt        
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}	
                % Use the ladder structure (much more efficient)
                xhi = dfbrec_l(Watermarked_image{end}, dfilt);
                
            otherwise	
                % General case
                xhi = dfbrec(Watermarked_image{end}, dfilt); 
        end
        
        x = lprec(xlo, xhi, h, g);
   
    else    
        % Special case: length(y{end}) == 3
        % Perform one-level 2-D critically sampled wavelet filter bank
        watermarked_image = wfb2rec(xlo, Watermarked_image{end}{1}, Watermarked_image{end}{2}, Watermarked_image{end}{3}, h, g);
    end
end
D=cover_object-watermarked_image;
%weight correction
alpha=0.5;
S=D*alpha;
watermarked_image1=watermarked_image+S;
watermarked_image_int=uint8(watermarked_image1);
path='C:\Users\USER\Desktop\';
iframe=1;
fname=strcat(path,int2str(iframe),'.jpg');
disp(fname);
imwrite(watermarked_image_int,fname,'jpg');
figure(1);
imshow(watermarked_image_int,[]);
title('Watermarked Image');
