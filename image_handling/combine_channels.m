function [combined_images] = combine_channels(images,noise)
%multicoil reconstruction
%   Detailed explanation goes here
% 
% noise = reshape(noise,numel(noise)/2, 2);
% noise = permute(noise,[2 1]);
% M = size(noise,2);
% Rn = (1/(M-1))*(noise*noise');
% Bn = (mean(abs(noise(:)))^2)./abs(noise(length(noise)./2))^2;
% Rnscaled = (Rn);


eta = noise;

n_channels = size(images,6);
psi = [n_channels,n_channels];
dims = size(images);
channel_data = (reshape(images,[],n_channels));
psi = (1/(length(noise)-1))*(eta*eta');
L = chol(psi,'lower');
L_inv =inv(L);
%  L_inv =1;
data_scaled = (L_inv*permute(channel_data,[2 1]));
channel_data = reshape(data_scaled',dims);
% channel_data(:,:,:,1) = 0;
% channel_data(:,:,:,2) = 0;
% channel_data(:,:,:,3) = 0;
% channel_data(:,:,:,4) = 0;
% channel_data(:,:,:,5) = 0;
% channel_data(:,:,:,6) = 0;
combined_images = sqrt(sum(channel_data.*conj(channel_data),6));


% for n=1:size(images,6)
% coil_images(:,:,:,n) = coil_images(:,:,:,n).*(noise(1)./noise(n));
% end
% 
% combined_images = sqrt(sum(coil_images.*conj(coil_images),4));