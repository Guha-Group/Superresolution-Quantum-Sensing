function ax = imagescComplex(xlim,ylim,field,colorscale,colorwheel_flag)
% Plots an image of complex matrix where the color (hue) represents the
% phase and the brightness represents the amplitude. Adapted from phplot by
% Iacopo Mochi.
%%%%%%%%%%%%%%%%%%%%
%%%%%% INPUTS %%%%%%
%%%%%%%%%%%%%%%%%%%%
%   xlim            : [1,2] vector of x-coordinate limits for image
%   ylim            : [1,2] vector of y-coordinate limits for image
%   field           : MxN complex matrix
%   colorscale        : parameter between (0,1] for changing phase coloring 
%                     saturation. A lower value saturates the color more.
%   colorwheel_flag : boolean flag for showing color wheel
%%%%%%%%%%%%%%%%%%%%%
%%%%%% INPUTS %%%%%%
%%%%%%%%%%%%%%%%%%%%
% Author: Nico Deshler
Im=imag(field);
Re=real(field);
phase=atan2(Im,Re);
amplitude=abs(field);
amplitude=amplitude/max(amplitude(:));
amplitude=amplitude.^(colorscale); % limit/boost color saturation via exponentiation

A=zeros(size(field,1),size(field,2),3);     %Declare RGB array
A(:,:,1)=0.5*(sin(phase)+1).*amplitude;     %Red
A(:,:,2)=0.5*(sin(phase+pi/2)+1).*amplitude;%Green
A(:,:,3)=0.5*(-sin(phase)+1).*amplitude;    %Blue
ax = imagesc(xlim,ylim,A);
axis square

if colorwheel_flag
    [X,Y] = meshgrid(linspace(-1,1,1001));
    [T,R] = cart2pol(X,Y);
    R(R > 1) = nan; 
    Z = R.*exp(1i*T);
    figure
    imagescComplex([-.5,.5],[-.5,.5],Z,1,0);
    axis square; axis off;
end

%{
if colorbar_flag
    s=-pi:2*pi/255:pi;
    c(:,1)=0.5*(sin(s)+1);     %Red
    c(:,2)=0.5*(sin(s+pi/2)+1);%Green
    c(:,3)=0.5*(-sin(s)+1);    %Blue
    colormap(c);
    cbar = colorbar;
    cbar.Limits = [0,1];
    cbar.Ticks = [0,.5,1];
    cbar.TickLabels={'$-\pi$','$0$','$\pi$'};
    cbar.TickLabelInterpreter='latex';
end
%}

end