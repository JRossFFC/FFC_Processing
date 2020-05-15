function [windowed_kspace] = windowkspace(kspace,window_size,windowtype)
%Windows kspace with a window_size point window with function windowtype
%   Detailed explanation goes here

dims = size(kspace); %we can reshape back to this dimensionality later
kspace = reshape(kspace, dims(1), dims(2),[]); %reshape to a 3D matrix for easier processing

windowdim1 = round(window_size*dims(1)); %ensure the window is integer sized
if mod(windowdim1,2)
    windowdim1 = windowdim1 + 1; %ensure window is symmetric
end
windowdim2 = round(window_size*dims(2));
if mod(windowdim2,2)
    windowdim2 = windowdim2 + 1;
end

switch windowtype %handle the string for user choice of windowing function
    case 'None'
        windowtype = 'rectwin';
        opts = [];
    case 'Hann'
        windowtype = 'hann';
        opts = [];
    case 'Blackman'
        windowtype = 'blackman';
        opts = [];
    case 'Hamming'
        windowtype = 'hamming';
        opts = [];
    case 'Kaiser-Bessel'
        windowtype = 'kaiser';
        opts = 3;
    case 'Tukey'
        windowtype = 'tukeywin';
        opts = 0.5;
    otherwise
        windowtype = 'rectwin';
        opts = [];
end

if isempty(opts)    %deal with windows that need additional inputs eg kaiser
    window1 = window(windowtype,windowdim1);
    window2 = window(windowtype,windowdim2);
else
    window1 = window(windowtype,windowdim1,opts); 
    window2 = window(windowtype,windowdim2,opts);
end

windowfull = window1*window2';  %create a 2D window

%deal with cases where the window is set to be larger or smaller than the kspace matrix
if length(windowfull)>size(kspace,1)    
    windowfull =  windowfull(round((size(windowfull,1) - size(kspace,1))/2)+1:end - round((size(windowfull,1) - size(kspace,1))/2),round((size(windowfull,2) - size(kspace,2))/2)+1:end - round((size(windowfull,2) - size(kspace,2))/2));
end

if length(windowfull)<size(kspace,1)
    windowfull = padarray(windowfull,[fix((size(kspace,1)-size(windowfull,1))/2), fix((size(kspace,2)-size(windowfull,2))/2)]);
end

if length(windowfull)<size(kspace,1)
   
    windowfull = padarray(windowfull,[1, 1],'post');
end

windowed_kspace = reshape(kspace.*windowfull,dims); %reshape back to the original matrix size

end

