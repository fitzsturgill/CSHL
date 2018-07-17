
function [X,Y,intercept_true,beta_true,RMSE,L1,auc_struct,line_handles ] = lasso_coord_desc_test( input_params )

%{

lasso_coord_desc_test
    Runs multivariate, 



Status

    - seems to work reasonably well with test dataset.
    - on NO intercept_hat case, doesn't seem to meet all the kkt criteria.
    - on intercept_hat case, it does seem to meet the criteria

    TODO :: 
        - Continue to fix multi-sample case.
            - Fix ktt tests - scale by n_samples?
            - Fix intercept_hat case - beta0 not properly defined.
            - Fix RMSE calculations below.
        - Figure out how to do cross-validation.

%}

%% Default parameters

% Standard parameters
default_params = struct();
default_params.min_lambda     = 0.001;
default_params.max_lambda     = 5;
default_params.n_lambdas      = 20;
default_params.lambda_linlog  = 'log'; % 'Linear' or 'log' defines spacing of lambda.
default_params.n_predictors   = 10;
default_params.n_dimensions   = 21;
default_params.n_samples      = 50;
default_params.SNR            = 1;
default_params.standardize    = 0;    % Rotate basis vectors?
default_params.regularization = 2;    % Acts on every second ß :: { 0 [do nothing] , 1 [shrink by 5x], 2 [set to zero] }
default_params.do_rotation    = 0;    % Rotate basis vectors?
% Over-ride internal variables - these are generated below if not passed.
default_params.X              = NaN;  % Placeholder for "randomize me please".
default_params.Y              = NaN;  % Placeholder for "randomize me please".
default_params.intercept_true = NaN;  % Placeholder for "randomize me please".
default_params.beta_true      = NaN;  % Placeholder for "randomize me please".
% AUC edge detector
% 'double_baseline_error' :: Samples the range from error(OLS) to 2x error(OLS)
% 'midpoint_RMSE'         :: Samples min(RMSE) to mean( [ min(RMSE), max(RMSE) ])
% 'midpoint_L1'           :: Samples min(RMSE) to RMSE @ mean( [ min(L1), max(L1) ])
% 'peak_inflection_point' :: Finds inflection by 2nd derivative, takes biggest peak
% 'last_inflection_point' :: Finds inflection by 2nd derivative, takes last peak
% 'cross-validation'      :: NOT IMPLEMENTED :: Find lambda with highest cross-val score.
default_params.auc_edge_detector = 'midpoint_RMSE';
% Can pass the edge of the AUC if defined from previous runs.
default_params.auc_edge_RMSE  = NaN;  % Placeholder for "calculate me please".
default_params.auc_edge_L1    = NaN;  % Placeholder for "calculate me please".
% Misc
default_params.do_figures     = 1;
default_params.reset_figures  = 1;
default_params.fig_suffix     = '';  % Can pass in a suffix for the figure titles.

default_params.verbose        = 1;
default_params.solver = 'lasso_coord_desc_test'; %By definition.

% Parse parameters and place in main namespace
if nargin > 0
    parsed_params = parse_defaults( default_params, input_params);
else 
    parsed_params = default_params;
end
fn = fieldnames(parsed_params);
for i = 1:length(fn)
    eval(sprintf('%s = parsed_params.%s;',fn{i},fn{i}));
end

%% Set defaults.

if ~exist('do_rotation','var')
    do_rotation = 0;
end

% beta_true is a column vector
% If passed, don't regularize it further.
% If X and Y are passed variables, leave beta_true as a NaN.
if all( [any(isnan(beta_true)) any(isnan(X)) any(isnan(Y))])
    beta_true = (rand(n_predictors,1)-0.5)*5;
    % Shrink some beta_true, or set some to 0, or do nothing.
    switch regularization
        case 0
            %beta_true = beta_true;
        case 1
            beta_true(1:2:end) = beta_true(1:2:end)/5;
        case 2
            beta_true(1:2:end) = 0;
    end
end

% Include an intercept term?  If intercept is NaN, we don't attempt.  If
% any other value, we attempt to find it.
% beta_plot is used for plotting.
if isnan(intercept_true) 
    use_intercept = 0;
    beta_plot = beta_true;
else
    error('Intercept not supported yet. Ask AGV to fix it.')
    use_intercept = 1;
    beta_plot     = [intercept_true; beta_true];
end

% X :: Basis vectors 
% We allow this to be passed in so that it can stay consistent.
if isscalar(X)
    assert(isnan(X),'X matrix is not properly formed - should either be a full matrix of predictors or a flag variable NaN for randomization.')
    X = rand(n_dimensions,n_predictors);
else
    n_predictors = size(X,2);
    n_dimensions = size(X,1);
end

% Z-score predictors?
if standardize
    % Allow nans in predictor matrix - they will be set to 0 after
    % standardization.
    X = nanzscore(X); 
    X(isnan(X)) = 0;
    Y = bsxfun(@minus,Y,mean(Y,1));
    Y(isnan(Y)) = 0;
end

% Generate Y from original basis vectors with noise.
% We allow this to be passed in so that it can stay consistent.
if isscalar(Y) && isnan(Y)
    if use_intercept
        Y = bsxfun(@plus, X*beta_true + intercept_true , (rand(n_dimensions,n_samples)-0.5)/(SNR/sum(beta_true)) );
    else
        Y = bsxfun(@plus, X*beta_true, (rand(n_dimensions,n_samples)-0.5)/(SNR/sum(beta_true)) );
    end
else
    n_samples = size(Y,2);
end

% A (column) rotation of X is a rotation of the basis vectors
if do_rotation
    X = random_rotation(X);
end

% Define range of lambda and include 0 if not otherwise present.
switch lambda_linlog
    case 'log'
        lambdas = logspace(log10(max(min_lambda,0.001)),log10(max_lambda),n_lambdas);
    case 'linear'
        lambdas = linspace(min_lambda,max_lambda,n_lambdas);
    otherwise
        error('I don''t recognize the lambda_linlog (%s) as an option.',lambda_linlog)
end

%Ensure that min(lambdas) == 0 to give the OLS estimator.
if min(lambdas)~=0,
    lambdas = [0 lambdas]; 
    n_lambdas = length(lambdas); 
end

if not(do_figures)
    line_handles = [];
end

% Structures to store results of LASSO
lasso_results = cell(n_lambdas,1);
beta_hat_matrix = nans(n_lambdas,n_predictors+use_intercept);


%% Iterate over lasso candidates and calculate results.

for i = 1:n_lambdas
    
    if verbose,
        fprintf('Iteration %3.0f :: lambda %.3f :: ',i,lambdas(i))    
    end

    lasso_options.lambda    = lambdas(i);
    lasso_options.maxiter   = 200;
    lasso_options.use_intercept = use_intercept;
    lasso_options.verbose   = verbose;    
    %lasso_options.tol       = 1E-9;

    lasso_results{i}        = lasso_coord_desc( X, Y, lasso_options);
    
    if use_intercept
        beta_hat_matrix(i,1) = lasso_results{i}.intercept_hat;
        beta_hat_matrix(i,2:end) = lasso_results{i}.beta_hat;
    else
        beta_hat_matrix(i,:) = lasso_results{i}.beta_hat;
    end
end

% Calculate mean RMSE across samples, as well as 1-norm of ß coefficients.
if use_intercept
    RMSE = sqrt(  mean(  bsxfun(@minus,Y,bsxfun(@plus,X*beta_hat_matrix(:,2:end)',beta_hat_matrix(:,1)')).^2));    
else

    % Original code for ONE__SAMPLE case.
    % RMSE = sqrt(mean( bsxfun(@minus,Y,X * lr').^2));
    % Y     :: [ n_dimensions x n_predictors  ]  
    % X     :: [ n_dimensions x n_predictors  ]
    % lr    :: [ n_lambdas    x  n_predictors ] 
    % X*lr' :: [ n_dimensions x n_lambdas     ] 
    [RMSE, RMSE_std]= deal(nans(1,n_lambdas));
    for i = 1:(n_lambdas)
        % Mean RMSE across samples.
        % For each lambda:
        % 1) Calculate squared error
        % 2) Take mean across dimensions
        % 3) Take sqrt
        % 4) Take mean across samples.
        RMSE(i)     =  mean(sqrt(mean(bsxfun(@minus,Y,X * beta_hat_matrix(i,:)').^2)));
        assert(RMSE(i) >= 0,'RMSE is somewho magically negative');
        
        RMSE_std(i) =  std( sqrt(mean(bsxfun(@minus,Y,X * beta_hat_matrix(i,:)').^2)));
    end
end

% Calculate L1 of coefficients.
L1 = sum(abs(beta_hat_matrix'));

%% Plot 1 :: Regression as a function of lambda

if do_figures
    
% Useful handles to pass back to calling function.
line_handles = {};

lambdas_to_plot = lambdas;
lambdas_to_plot(lambdas == 0) = min(lambdas(2:end))*0.8;

namedFigure(sprintf('Regression analysis :: %s',fig_suffix))
if reset_figures, clf, end
%namedFigure(sprintf('Regression analysis :: do_rotation = %.0f',do_rotation))
colormap red_white_blue

% Common color axis.
beta_axis = max(abs(beta_hat_matrix(:))) * [-1 1] + [-1 1];
colormap red_white_blue
% Padding for pcolor plots
pcolor_lambdas    = [lambdas lambdas(end)];
pcolor_predictors = (1:n_predictors+1)-0.5;
[meshx,meshy] = meshgrid( pcolor_lambdas, pcolor_predictors);
[pcolor_lr, pcolor_error] = deal(zeros( size(beta_hat_matrix') + [1 1]));
pcolor_lr(1:size(beta_hat_matrix,2),1:size(beta_hat_matrix,1)) = beta_hat_matrix';
pcolor_error(1:size(beta_hat_matrix,2),1:size(beta_hat_matrix,1)) = bsxfun(@minus,beta_hat_matrix,beta_plot')';

%  Plot regressio coefficients as heatmap.
if ~do_rotation % Don't over-write if not rotating.
    if ~isnan(beta_true)
        subplot(10,10,[1:9 11:19])
    else
        splot = bsxfun(@plus,[1:9],[0:10:30]');
        subplot(10,10,sort(splot(:)))
    end
    pcolor(meshx,meshy,pcolor_lr); shading flat;
    set(gca,'XScale','log')
    set(gca,'XTick',[])
    axis ij
    caxis(beta_axis)
    colorbar off
    title('Estimated')
    ylabel('Predictors')
end

% Plot beta_true from original model (if known).
if ~isnan(beta_true)
    subplot(10,10,[10 20])
    imagesc(beta_plot)
    caxis([-max(abs(beta_plot))-1 max(abs(beta_plot))+1])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    title('Actual')
    axis ij
    colorbar
    axis off
    
    % Error - estimated beta_true - real beta_true;
    subplot(10,10,[21:29 31:39])
    pcolor(meshx,meshy, pcolor_error); shading flat;
    set(gca,'XScale','log')
    set(gca,'XTick',[])
    caxis(beta_axis)
    ylabel('Predictors')
    text(lambdas(2)*1.1,2,'Error : beta - beta_hat','Interpreter','none')
    axis ij
end

% Estimated beta_true coefficients.
subplot(10,10,[41:49 51:59 61:69])
hold on
if do_rotation == 0
    line_handles{end+1} = plot(lambdas_to_plot,beta_hat_matrix,'-','LineWidth',2);
    line_handles{end+1} = plot(lambdas_to_plot(1),beta_plot,'o','MarkerSize',10,'LineWidth',2);;
else
    line_handles{end+1} = plot(lambdas_to_plot,beta_hat_matrix,':','LineWidth',1);
end
%hline(beta_plot,'k:')
%plot(lambdas_to_plot(1),beta_plot,'-.')
%legend( cellfun( @(x) num2str(x), num2cell(1:length(beta_plot)), 'UniformOutput',false),'Location','EastOutside')
ylim(beta_axis)
xlim([min(lambdas_to_plot) max(lambdas_to_plot)] .* [0.8 1])
set(gca,'XScale','log')
set(gca,'XTick',[])
box on

% RMSE and L1-norm of coefficients.
subplot(10,10,[71:79 81:89 91:99])
hold on
[hAx,hLine1,hLine2] = plotyy(lambdas_to_plot,RMSE,...
                             lambdas_to_plot,sum(abs(beta_hat_matrix')),'semilogx','semilogx');
line_handles{end+1} = hLine1;
line_handles{end+1} = hLine2;
xlabel('lambda')
set(hLine1,'LineStyle','-','LineWidth',2)
%set(hLine2,'LineStyle','o')
h_error = errorbar(lambdas_to_plot,RMSE,RMSE_std);
set(gca,'XScale','log')
set(hAx(1),'Xlim',([min(lambdas_to_plot) max(lambdas_to_plot)] .* [0.8 1]))
set(hAx(2),'Xlim',([min(lambdas_to_plot) max(lambdas_to_plot)] .* [0.8 1]))
%set(hAx(1),'Ylim',[0 2]);
%set(hAx(2),'Ylim',[0 12]);
set(hAx(2),'XScale','log')
ylabel(hAx(1),({'RMSE (mean +/- stdev)',sprintf('for %.0f samples',n_samples)}))
ylabel(hAx(2),{'1-Norm of','beta_true Coefficients'})
%linkaxes(hAx)
%axis tight
%ylims = ylim;
%ylim([0 max(ylims(2),1.1)]);

end
%% Plot 2 :: RMSE vs. L1

if do_figures
    
namedFigure(sprintf('RMSE vs. L1 :: %s',fig_suffix))
if reset_figures, clf, end


logscale_fig2 = 1;
if logscale_fig2
    set(gca,'XScale','log')
    set(gca,'YScale','log')
end

hold on

if do_rotation == 0
    
    line_handles{end+1} = plot(RMSE,L1,'k--o','LineWidth',2);
    
    
    % Log
    if logscale_fig2
        % Log
        xmin = min(RMSE)*.66;
        xmax = max(RMSE)*1.5;
         %ylim([ylim .* [1 1.25] + [max([1E-4,min(L1)]) 0] ])
    else
        % Linear
        ylim([ylim .* [1 1.25]])
        xmin =  max([0,min(RMSE)-0.66*range(RMSE)]);
        xmax =  max(RMSE) + max([min([ 1.5 * range(RMSE), min(RMSE)-xmin]),range(RMSE)/10]);
    end
    xlim([ xmin xmax]);
    ylims = ylim;
    xlims = xlim;
    
    vline(RMSE(1),'k:')
    vline(RMSE(end),'k:')
    hline(L1(1),'k:')
    
    % Plot lines on the right side.
    for label = [0, 0.01, 0.1, 1]
        
        plot( [ xlims(end)-0.05*diff(xlims) xlims(end)], ones(1,2) .* L1(find(lambdas>=label,1,'first')) ,'k','LineWidth',2);
        text( xlims(2)*0.95, L1(find(lambdas>=label,1,'first'))+diff(ylims)/50  ,sprintf('\\lambda = % 1.1g',label))
    end

    plot(RMSE(1),L1(1),'ko','MarkerSize',10)
    
else
    hold on
    plot( RMSE,    L1,   'r--.','LineWidth',  1)
    plot( RMSE(1), L1(1),'ro',  'MarkerSize',10)
end

%hline(L1(find(lambdas>1,1)),'g:')
xlabel('RMSE')
ylabel('L1')
axis on

end


%% Calculate AUC 

% If AUC edge is not defined (ie., passed as an input parameter) it will be
% a NaN.  In that case, calculate it.
if isnan(auc_edge_RMSE) || isnan(auc_edge_L1) 
    
    % We want to calculate the AUC of the RMSE-L1 curve, but we know that
    % the right side of this curve is dominated by the size of the last
    % remaining predictors and thus want to ignore it.  There are several
    % methods for defining where to place the right hand side of this
    % curve.
    %
    % 'full_range'            :: Samples the full from Lambda 0 .. max.
    % 'double_baseline_error' :: Samples the range from error(OLS) to 2x error(OLS)
    % 'midpoint_RMSE'         :: Samples min(RMSE) to mean( [ min(RMSE), max(RMSE) ])
    % 'midpoint_L1'           :: Samples min(RMSE) to RMSE @ mean( [ min(L1), max(L1) ])
    % 'peak_inflection_point' :: Finds inflection by 2nd derivative, takes biggest peak
    % 'last_inflection_point' :: Finds inflection by 2nd derivative, takes last peak
    % 'first_point'           :: NOT IMPLEMENTED :: Samples the first point of beta - ie., the L1 norm of ß at lambda = 0.
    % 'cross-validation'      :: NOT IMPLEMENTED :: Find lambda with highest cross-val score.
    %
    % We define two variables: 
    %  AUC_edge_ind   :: The last index of RMSE/L1 vectors that's in the AUC box.
    %  RMSE_edge_trapz :: The actual edge, which may lie between values of the RMSE/L1 vectors.  
    %  L1_edge_trapz ::   The L1 value at the last point. This is calculate by spline interpolation
    %                     and included as the last point in the AUC calculation by trapz()

    
    spline_inds = 1:find(diff(RMSE)==0,1); % Used in spline estimates, since we can't repeat RMSE values.
    switch auc_edge_detector
        
      
        case 'full_range' 
            auc_edge_ind  = max(spline_inds);
            auc_edge_RMSE = RMSE(auc_edge_ind);
            auc_edge_L1   = L1(auc_edge_ind);
        
        case 'double_baseline_error'
            % Take values up to 2x the OSL RMSE (which is included in the RMSE vector)
            auc_edge_ind  = find(RMSE < min(RMSE)*2, 1,'last');
            auc_edge_RMSE = min(RMSE)*2;
            auc_edge_L1   = spline( RMSE(spline_inds), L1(spline_inds), auc_edge_RMSE);

        case 'midpoint_RMSE'
            % Limit our AUC box to 50% of the range of RMSE observed for lambda {0...max}
            auc_edge_ind  = find(RMSE < mean([min(RMSE),max(RMSE)]),1,'last');
            auc_edge_RMSE = mean([min(RMSE),max(RMSE)]);
            auc_edge_L1   = spline( RMSE(spline_inds), L1(spline_inds), auc_edge_RMSE);
            
        case 'midpoint_L1'
            % Limit our AUC boxx to 50% of the range of observed L1 error.
            % Note that L1 ~ 1/lambda so is descending.
            auc_edge_ind  = find(L1 > mean([min(L1),max(L1)]),1,'last');
            auc_edge_L1   = mean([min(L1),max(L1)]);
            auc_edge_RMSE = spline( L1(spline_inds), RMSE(spline_inds), auc_edge_L1);
            
        case {'peak_inflection_point','last_inflection_point'}
            % Generate splines for ~high-resolution peak finding. 
            RMSE_spline = linspace(RMSE(1),RMSE(spline_inds(end)),100);           % Generate evenly spaced RMSE vector.
            L1_spline = spline( RMSE(spline_inds), L1(spline_inds), RMSE_spline); % Interpolate L1 from spline..
            % Find peaks in the second derivative of spline-interpolated L1.
            [peaks, peak_locs] = findpeaks(  diff(diff(L1_spline)) );
            plot(RMSE_spline(peak_locs+2),L1_spline(peak_locs+2),'ro') %Plot all peaks.
            switch auc_edge_detector
                case 'peak_inflection_point'
                    auc_edge_RMSE = RMSE_spline(peak_locs(peaks == max(peaks))+2); % Find RMSE at best peak.
                    auc_edge_L1   = L1_spline(  peak_locs(peaks == max(peaks))+2); % Find RMSE at best peak.
                    auc_edge_ind  = find(RMSE < auc_edge_RMSE,1,'last');
                case 'last_inflection_point'
                    auc_edge_RMSE = RMSE_spline(peak_locs(end));
                    auc_edge_L1   = L1_spline(  peak_locs(end));
                    auc_edge_ind  = find(RMSE < auc_edge_RMSE,1,'last');
            end
        otherwise
            error('Method for auc_edge_detector not found :: %s',auc_edge_detector)
    end

    hline( auc_edge_L1,   'r:')
    vline( auc_edge_RMSE, 'r:')
    vline( min(RMSE),     'k:')
    vline( max(RMSE),     'k:')
    
else
    % Both auc_edge_RMSE and auc_edge_L1 have been passed, so we assume
    % these.
    % Calculate new index for the right-most sampled point in the AUC.
    % We will use auc_edge_RMSE and auc_edge_L1 as passed parameters;
    auc_edge_ind  = find(RMSE < auc_edge_RMSE,1,'last');
end
   
% Calculate the AUC for the RMSE/L1 curve up to the calculated edge.
auc = trapz( [ RMSE(1:auc_edge_ind)  auc_edge_RMSE ] , ...
             [   L1(1:auc_edge_ind)    auc_edge_L1 ] );

%% Output variables

auc_struct = struct();
auc_struct.auc_edge_RMSE = auc_edge_RMSE;
auc_struct.auc_edge_L1 = auc_edge_L1;
auc_struct.auc = auc;

drawnow

%     %%
%     figure; 
%     subplot(1,2,1); hold on;
%     plot(RMSE_spline(3:end),diff(diff(L1_spline)))
%     plot(RMSE_spline(peak_locs+2),peaks,'ro')
%     subplot(1,2,2); hold on
%     plot(RMSE_spline(3:end),L1_spline(3:end))
%     plot(RMSE_spline(peak_locs+2),L1_spline(peak_locs+2),'ro')
%     %%
%     keyboard
%     
    

