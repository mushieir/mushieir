clc;
clear all;
close all;
warning off all;
tic;

%% �������飬����ͼ����������
I0=zeros(512,512,512);
temp=zeros(512,512);
%����������ɨ% 
PicDic=1:1023;
for i=1:512
    %���������������
    temp=double(imread(strcat('E:\polarization program\201804_05experiment\ex0329\New folder\20180327_5_Camera_',num2str(PicDic(1024-i)+2000,'%04d'),'.tif')));
    for j=1:i
        I0(:,512-i+j,j)=temp(:,j);    
    end
end
for i=513:1013
    %���������������
     temp=double(imread(strcat('E:\polarization program\201804_05experiment\ex0329\New folder\20180327_5_Camera_',num2str(PicDic(1024-i)+2000,'%04d'),'.tif'))); 
    for j=i-511:512    
        I0(:,512-i+j,j)=temp(:,j);
    end
end

for y=1:1:256;
    I_flip(y+256,:,:)=I0(257-y,: ,:);
    I_flip(y,:,:)=I0(y+256,:,:);
end
I0=I_flip;
% save('I0_0412_0.mat','I0');

%% ȥ����
I0_bk=zeros(512,512,512);
load YFilter
for i=1:512
    S=fftshift(fft2(reshape(I0(i,:,:),512,512)));                                     %Ƶ��任
    I0_bk(i,:,:)=real(ifft2(fftshift(S.*YFilter)));
end
save('I0_bk_0412_0.mat','I0_bk');


%% ƫ��������������帴ԭ��ȡ
S0=zeros(512,512,137);
S1=zeros(512,512,137);
S2=zeros(512,512,137);
S3=zeros(512,512,137);
temp2=zeros(1,1,512);
temp4=zeros(1,1,512);
temp6=zeros(1,1,512);
for ii=1:512
    for jj=1:512
        temp2(1,1,219:292)=I0_bk(ii,jj,219:292); 
        temp3=2*abs(ifft(fftshift(temp2)));
        S0(ii,jj,:)=temp3(120:256);     
        temp4(1,1,366:438)=I0_bk(ii,jj,366:438);
        temp5=2*abs(ifft(fftshift(reshape(temp4,1,512))));   % no phase corr 
        S1(ii,jj,:)=temp5(120:256);       
        temp6(1,1,293:365)=I0_bk(ii,jj,293:365);
        temp7=abs(ifft(fftshift(reshape(temp6,1,512))));   % no phase corr 
        S2(ii,jj,:)=temp7(120:256);
        S3(ii,jj,:)=0.25*temp7(120:256);   % no phase corr 
    end
end
save('S0_test.mat','S0');
save('S1_test.mat','S1');
save('S2_test.mat','S2');
save('S3_0330_1.mat','S3');

%% ���ɫ��ԭ
load CIE1931rgb       %������������飬�ֱ�ΪCIEr��CIEg��CIEb,ÿ����Ϊ1*81��double��ÿ��������380nm��ʼ���5����780nm���ף�
% ��������137�������1/480����һ������1/980,ʵ��1/480�����ڵ�m�㣬���±ߵ�137�滻��m����

S0_Cim=zeros(512,512,3);
S1_Cim=zeros(512,512,3);
S2_Cim=zeros(512,512,3);
S3_Cim=zeros(512,512,3);
wavenum=1/960:1/960/(137-1):1/480;  
wavelength=1./wavenum;
AA=find(wavelength<=780);
XX=137:-1:AA(1);
number=length(XX);
wavelength(wavelength>780)=[];
CIEWaveLength=380:5:780;
CIER=interp1(CIEWaveLength,CIEr,wavelength,'spline');
CIEG=interp1(CIEWaveLength,CIEg,wavelength,'spline');
CIEB=interp1(CIEWaveLength,CIEb,wavelength,'spline');
for ii=1:512
    for jj=1:512
        S0_Cim(ii,jj,1)=sum(reshape(S0(ii,jj,XX),number,1).*CIER(end:-1:1)');
        S0_Cim(ii,jj,2)=sum(reshape(S0(ii,jj,XX),number,1).*CIEG(end:-1:1)');
        S0_Cim(ii,jj,3)=sum(reshape(S0(ii,jj,XX),number,1).*CIEB(end:-1:1)');
        
        S1_Cim(ii,jj,1)=sum(reshape(S1(ii,jj,XX),number,1).*CIER(end:-1:1)');
        S1_Cim(ii,jj,2)=sum(reshape(S1(ii,jj,XX),number,1).*CIEG(end:-1:1)');
        S1_Cim(ii,jj,3)=sum(reshape(S1(ii,jj,XX),number,1).*CIEB(end:-1:1)');
        
        S2_Cim(ii,jj,1)=sum(reshape(S2(ii,jj,XX),number,1).*CIER(end:-1:1)');
        S2_Cim(ii,jj,2)=sum(reshape(S2(ii,jj,XX),number,1).*CIEG(end:-1:1)');
        S2_Cim(ii,jj,3)=sum(reshape(S2(ii,jj,XX),number,1).*CIEB(end:-1:1)');
        
        S3_Cim(ii,jj,1)=sum(reshape(S3(ii,jj,XX),number,1).*CIER(end:-1:1)');
        S3_Cim(ii,jj,2)=sum(reshape(S3(ii,jj,XX),number,1).*CIEG(end:-1:1)');
        S3_Cim(ii,jj,3)=sum(reshape(S3(ii,jj,XX),number,1).*CIEB(end:-1:1)');
    end
end
S0_Cim=S0_Cim./max(max(max(S0_Cim)));
S1_Cim=S1_Cim./max(max(max(S1_Cim)));
S2_Cim=S2_Cim./max(max(max(S2_Cim)));
S3_Cim=S3_Cim./max(max(max(S3_Cim)));
save('S0_Cim.mat','S0_Cim');
save('S1_Cim.mat','S1_Cim');
save('S2_Cim.mat','S2_Cim');
save('S3_Cim.mat','S3_Cim');

%% ͼ����ǿ
S0_T=rgb2hsv(S0_Cim);
S1_T=rgb2hsv(S1_Cim);
S2_T=rgb2hsv(S2_Cim);
S3_T=rgb2hsv(S3_Cim);      %�任��hsv�ռ�
T0=zeros(512,512,3);
T1=zeros(512,512,3);
T2=zeros(512,512,3);
T3=zeros(512,512,3);
L=ones(512,512);
L(1:230,256:257)=0;
L(283:512,256:257)=0;                  %�����ݲ��˲���
for i=1:3
    T0(:,:,i)=real(ifft2(fftshift(fftshift(fft2(S0_T(:,:,i))).*L)));
    T1(:,:,i)=real(ifft2(fftshift(fftshift(fft2(S1_T(:,:,i))).*L)));
    T2(:,:,i)=real(ifft2(fftshift(fftshift(fft2(S2_T(:,:,i))).*L)));
    T3(:,:,i)=real(ifft2(fftshift(fftshift(fft2(S3_T(:,:,i))).*L)));
end
T0(:,:,3)=I0(:,:,100)./max(max(I0(:,:,100))); %������ͼ������
S0_T=hsv2rgb(T0);
S1_T=hsv2rgb(T1);
S2_T=hsv2rgb(T2);
S3_T=hsv2rgb(T3);        %�任��rgb�ռ�

%% ƫ��ȣ�ƫ���ͼ��
n=21;
Aop=zeros(512,512);
Dop=zeros(512,512);
L1=ones(512,512);
L1(1:255,256:257)=0;
L1(256:512,256:257)=0;             %�����ݲ��˲���
Sum0=sum(S0(:,:,2:136),3);
Sum1=sum(S1(:,:,2:136),3);
Sum2=sum(S2(:,:,2:136),3);
Sum3=sum(S3(:,:,2:136),3);
for i=1:512
    for j=1:512
        Aop(i,j)=(1/2*atan(Sum2(i,j)./Sum1(i,j)))./pi*180;
        Dop(i,j)=sqrt(Sum1(i,j).^2+Sum2(i,j).^2)./Sum0(i,j);
    end
end
Dop=Dop/max(max(Dop));
% Aop=real(ifft2(fftshift(fftshift(fft2(Aop)).*L1)));
% Dop=real(ifft2(fftshift(fftshift(fft2(Dop)).*L1)));

%% �������
%%%%%%%%%%%%%%����Ҷ�ͼ%%%%%%%%%%%%%%%
figure(1)
imshow(I0(:,:,100),[]);
title('�Ҷ�ͼ')
%%%%%%%%%%%%%%���˹�п�˹��ɫͼ%%%%%%%%%%%%%
figure(2)
imshow(S0_T,[]);
title('S0')
figure(3)
imshow(S1_T,[]);
title('S1')
figure(4)
imshow(S2_T,[]);
title('S2')
figure(5)
imshow(S3_T,[]);
title('S3')
%%%%%%%%%%%%%���ƫ�����ƫ���%%%%%%%%%%%%%
figure(6)
imshow(Aop,[]);
title('ƫ���')
colormap hsv
colorbar
figure(7)
imshow(Dop,[]);
title('ƫ���')
%%%%%%%%%%%%�������������%%%%%%%%%%%%%%
figure(8)
I_c=I0(2:2:512,2:2:512,2:2:512);
[x,y,z]=meshgrid(2:2:512,2:2:512,2:2:512);             
xslice = [2,512]; yslice = [2,512]; zslice = [2,512];  %����ѡ��һЩ��Ƭ
GG=slice(x,y,z,I_c,xslice,yslice,zslice);    %permute��������
shading interp
axis off
title('�Ҷ�ͼ����������')
colormap jet     
colorbar
figure(9)
S0_c=S0(2:2:512,2:2:512,2:2:136);
[x,y,z]=meshgrid(2:2:512,2:2:512,2:2:136);             
xslice = [2,512]; yslice = [2,512]; zslice = [10,136];  %����ѡ��һЩ��Ƭ
GG=slice(x,y,z,S0_c,xslice,yslice,zslice);    %permute��������
shading interp
axis off
title('S0����������')
colormap jet     
colorbar
figure(10)
S1_c=S1(2:2:512,2:2:512,2:2:136);
[x,y,z]=meshgrid(2:2:512,2:2:512,2:2:136);             
xslice = [2,512]; yslice = [2,512]; zslice = [2,136];  %����ѡ��һЩ��Ƭ
GG=slice(x,y,z,S1_c,xslice,yslice,zslice);    %permute��������
shading interp
axis off
title('S1����������')
colormap jet     
colorbar
figure(11)
S2_c=S2(2:2:512,2:2:512,2:2:136);
[x,y,z]=meshgrid(2:2:512,2:2:512,2:2:136);             
xslice = [2,512]; yslice = [2,512]; zslice = [2,136];  %����ѡ��һЩ��Ƭ
GG=slice(x,y,z,S2_c,xslice,yslice,zslice);    %permute��������
shading interp
axis off
title('S2����������')
colormap jet     
colorbar
figure(12)
S3_c=S3(2:2:512,2:2:512,2:2:136);
[x,y,z]=meshgrid(2:2:512,2:2:512,2:2:136);             
xslice = [2,512]; yslice = [2,512]; zslice = [2,136];  %����ѡ��һЩ��Ƭ
GG=slice(x,y,z,S3_c,xslice,yslice,zslice);    %permute��������
shading interp
axis off
title('S3����������')
colormap jet     
colorbar

%%
toc;