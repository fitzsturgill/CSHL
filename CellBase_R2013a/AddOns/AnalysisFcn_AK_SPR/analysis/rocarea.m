function [D, P] = rocarea(x,y,varargin)
%
%  [D, P] = rocarea(x,y,{'boot',n},{'transform'},{'PLOT'})
%
%  Computes discriminability index or area under ROC curve. 
%
%  x, y : data
%  'boot', n: number of bootstraps
%  transform: 'swap' -- always gives you results between 0.5 - 1
%             'scale' -- scales to give you results from -1 to 1
%
% 'PLOT' -- plots the ROC curve
%
% AK 2/2002
% AK 4/2005

Nboot = 0;
TRANSFORM = 0;
if isempty(x) | isempty(y)
    D = NaN;
    P = NaN;
   return
end

if nargin > 2
    switch lower(varargin{1}) 
        case 'boot'
            Nboot = varargin{2};
        case 'swap'
              TRANSFORM = 1;
        case 'scale'
              TRANSFORM = 2;
        otherwise
              TRANSFORM = 0;
    end
    if nargin > 4
      switch lower(varargin{3}) 
        case 'swap'
              TRANSFORM = 1;
        case 'scale'
              TRANSFORM = 2;
      end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=x(:); y=y(:);


nbin = ceil(max(length(x)*1.2,length(y)*1.2)); % some automatic assignment for number of bins
MN = min([x;y]);    MX = max([x;y]);
bin_size = (MX-MN)/nbin;
bins = MN-bin_size:bin_size:MX+bin_size;
Lx = length(x); Ly = length(y);

%%%%%%%%%%%%%%%%%%%%
% ROC calculation
%%%%%%%%%%%%%%%%%%%%
D = auc(x,y,Lx,Ly,bins);


%%%%%%%%%%%%%%
% bootstrap
%%%%%%%%%%%%%
if Nboot > 0

 z = [x; y];
 for i=1:Nboot    
   order=round(rand(1,Lx+Ly)*(Lx+Ly-1))+1;  %resample
   px=z(order(1:Lx));           %resort
   py=z(order(Lx+1:Lx+Ly));
   Dboot(i)= auc(px,py,Lx,Ly,bins);     %recalculate D
 end

 %
 % Decide which side it should be on
 P = iprctile(Dboot,D);
 
 if D > mean(Dboot)
     P = 1 - P;
 end
else
   P = 1;
end


%%%%%%%%%%%%%%
% transform D
%%%%%%%%%%%%%
switch TRANSFORM
    case 1 
        D=abs(D-0.5)+0.5;       % 'swap'
    case 2
        D=2*(D-0.5);            % 'scale'
end

%%%%%%%%%%%%
%  P L O T
%%%%%%%%%%%
if nargin > 2 & strcmp(lower(varargin{end}), 'plot')
 figure(1)
  clf;hold on
  plot(cdf2,cdf1,'b','LineWidth',2);
  plot([0 1],[0 1],'k');
  xlabel('False alarm'); ylabel('Hit rate');
  title(num2str(D));
end

function D = auc(x,y,Lx,Ly,bins);

p = histc(x,bins);  q = histc(y,bins);

cdf1 = cumsum(p)/Lx;   cdf2 = cumsum(q)/Ly;

if isempty(cdf1) || isempty(cdf2)
    D = NaN;
else
    D=trapz(cdf1,cdf2);
end