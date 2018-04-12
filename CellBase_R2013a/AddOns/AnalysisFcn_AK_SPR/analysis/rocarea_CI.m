function varargout = rocarea_CI(x,y,varargin)
% 
%  [D, P, CI] = rocarea(x,y,{'boot',n},{'transform'},{'PLOT'})
%
%  Computes discriminability index or area under ROC curve. 
% D : discriminability index or auROC
% Confidence Interval : bootstrapped confidence interval for D statistic 
%  x, y : data
%  'boot', n: number of bootstraps
%  transform: 'swap' -- always gives you results between 0.5 - 1
%             'scale' -- scales to give you results from -1 to 1
%
% 'PLOT' -- plots the ROC curve
%
% AK 2/2002
% AK 4/2005
% FS 4/2018

Nboot = 0;
TRANSFORM = 0;
if isempty(x) || isempty(y)
    if nargout > 0
        varargout{1} = NaN; % D
    end
    if nargout > 1
        varargout{2} = NaN; % P
    end
    if nargout > 2
        varargout{3} = [NaN NaN]; % CI
    end    
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
% bootstrap- bootstrap P value and confidence interval for D statistic
%%%%%%%%%%%%%
if Nboot > 0

 z = [x; y]; % null hypothesis, x and y are drawn from same distribution
 DbootNull = NaN(Nboot, 1);
 Dboot = NaN(Nboot, 1);
 for i=1:Nboot    
   % can we reject null hypothesis?
   orderNull=round(rand(1,Lx+Ly)*(Lx+Ly-1))+1;  %resample null
   pxNull=z(orderNull(1:Lx));           %resort
   pyNull=z(orderNull(Lx+1:Lx+Ly));
   DbootNull(i)= auc(pxNull,pyNull,Lx,Ly,bins);     %recalculate D
   
   % confidence intervals for D
   orderx = round(rand(1,Lx)*(Lx-1))+1;  %resample x
   ordery = round(rand(1,Ly)*(Ly-1))+1;  %resample y
   px=x(orderx);           
   py=y(ordery);
   thisD = auc(px,py,Lx,Ly,bins);
   switch TRANSFORM
    case 1 
        thisD=abs(thisD-0.5)+0.5;       % 'swap'
    case 2
        thisD=2*(thisD-0.5);            % 'scale'
   end
   Dboot(i)= thisD;     % transformed
 end
    % P value for null hypothesis
    P = iprctile(DbootNull,D);
    if D > mean(DbootNull)
        P = 1 - P;
    end
    % compute confidence intervals from bootstrapped distributions
    CI = zeros(1,2);
    CI(1) = percentile(Dboot, .05);
    CI(2) = percentile(Dboot, .95);    
else
   P = 1;
   CI = [NaN NaN];
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
  p = histc(x,bins);  q = histc(y,bins);
  cdf1 = cumsum(p)/Lx;   cdf2 = cumsum(q)/Ly;
  hold on
  plot(cdf2,cdf1,'b','LineWidth',2);
  plot([0 1],[0 1],'k');
  xlabel('False alarm'); ylabel('Hit rate');
  title(num2str(D));
end
% set outputs
if nargout > 0
    varargout{1} = D;
end
if nargout > 1
    varargout{2} = P;
end
if nargout > 2
    varargout{3} = CI;
end
    

function D = auc(x,y,Lx,Ly,bins)

p = histc(x,bins);  q = histc(y,bins);

cdf1 = cumsum(p)/Lx;   cdf2 = cumsum(q)/Ly;

if isempty(cdf1) || isempty(cdf2)
    D = NaN;
else
    D=trapz(cdf1,cdf2);
end