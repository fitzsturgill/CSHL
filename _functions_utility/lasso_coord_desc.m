function beta_struct = lasso_coord_desc( X, Y, lasso_options)
% Multivariate lasso by coordinate descent, based off of the R
% implementation at
% http://jocelynchi.com/a-coordinate-descent-algorithm-for-the-lasso-problem/#the-multivariate-lasso-problem
%
%
%
% code_version = 0.06
%
% NOTES 
%
%   v.0.06
%       For multi-sample code, lambda seems way too strong. Scale by
%       n_samples?  Why is this arising?
%


% DEBUG PARAMETERS
% Allow break over iteration loop before maxiter?
allow_break = 1;

if nargin < 3,
    lasso_options = {};
end

% Default parameters
default_params = { ...
    'lambda',       0    , ...
    'tol',          1E-9 ,...
    'use_intercept',   0 ,...
    'maxiter',      5000 ,...
    'verbose',      true ...
    };
% Parse defaults.
parsed_params = parse_defaults(default_params,lasso_options);
fn = fieldnames(parsed_params);
for i = 1:length(fn)
    eval(sprintf('%s = parsed_params.%s;',fn{i},fn{i}));
end

% If no inputs, generate synthetic dataset.
switch nargin
    case 0
        % Generate synthetic dataset.
        fprintf('---------------------------------------------------\n')
        fprintf('Running lasso_coord_descent with synthetic dataset.\n')
        
        % Specify dimensions etc.
        generate_synthetic_dataset = true;
        n_dimensions  = 20;
        n_predictors  = 5;
        n_samples     = 6;
        intercept     = NaN; % Nan for "no intercept".  Otherwise will attempt.
        
        % 5 predictors in 10 dimensions.
        % Columns are independent predictors
        % Rows are dimensions.
        X = rand(n_dimensions,n_predictors);
        
        % beta is a column vector
        beta_true = (rand(n_predictors,1)-0.5)*5;
        beta_true(1:2:end) = 0;
        
        % Approximate SNR - this noise puts a ceiling on prediction accuracy.
        % You should probably compare to a standard OLS regression for sanity.
        SNR = 20;
        
        % Calculate respopnse
        % Only modeling additive noise right now.
        use_intercept = not(isnan(intercept));
        if use_intercept
            Y = bsxfun(@plus, X*beta_true + intercept , rand(n_dimensions,n_samples)/(SNR/sum(beta_true)) );
        else
            Y = bsxfun(@plus, X*beta_true , rand(n_dimensions,n_samples)/(SNR/sum(beta_true)) );
        end
    case 1
        error('Please specify both predictors and response variables.')
        
    otherwise
        % Utility variables
        generate_synthetic_dataset = false;
        n_dimensions = size(X,1);
        n_predictors = size(X,2);
        n_samples    = size(Y,2);
end

% Initialize estimates of beta etc.
beta_hat = zeros(size(X,2),1); % Beta is a column vector such that Y = X*beta.

% Sanity checks
% Ensure that predictors and outputs have the same number of dimensions.
assert(size(X,1) == size(Y,1)   ,'Mismatch in dimensions: Rows of X (dimensions) must match the length of Y.')
assert(size(X,2) == length(beta_hat),'Mismatch in number of predictors :: Columns of X (# predictors) must match the lenght of beta_hat.')
assert(~any(isnan(X(:))),'Nans found in predictor matrix X')
assert(~any(isnan(Y(:))),'Nans found in response matrix Y')


%%  Begin translation from R %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj         = nans(maxiter+1,1);
beta_list    = cell(maxiter+1,1);
beta_list{1} = beta_hat;
if use_intercept,
    
    %     intercept_list <- numeric(length(maxiter+1))
    %     intercept_hat <- sum(y-X%*%beta_hat)/(length(y))
    %     intercept_list[1] <- intercept_hat
    intercept_hat_list = nans(maxiter+1,1);
    % ONE__SAMPLE
    %intercept_hat = sum(Y - X * beta_hat) / n_dimensions;
    % MULTISAMPLE
    error('FIX intercept_hat HERE')
    intercept_hat = sum(mean(bsxfun(@minus,Y,X * beta_hat))) / n_dimensions;
    intercept_hat_list(1) = intercept_hat;
    
    for j = 1:maxiter
        
        for k = 1:n_predictors
            
            % r <- y - X[,-k]%*%beta_hat[-k] - intercept_hat*rep(1,n_dimensions)
            not_k = [1:k-1,k+1:n_predictors];
            
            % ONE__SAMPLE :: y - X*beta_hat is a vector of size [n_dimension,1]
            %                so r is size(y) b/c r = y - X(:,not_k) * beta_hat(not_k) - intercept_hat * ones(n_dimensions,1);
            % MULTISAMPLE :: We want the prediction error for each sample separately.
            %                so r = size(Y);
            r = bsxfun(@minus,Y,X(:,not_k) * beta_hat(not_k) - intercept_hat * ones(n_dimensions,1));
            
            % beta_hat[k] <- (1/norm(as.matrix(X[,k]),"F")^2)
            %            *   soft_thresholding(t(r)%*%X[,k],n_dimensions*lambda)
            % ONE__SAMPLE :: r' * X(:,k) is size [1,1]           <-- [1,n_dimensions] * [n_dimensions 1]
            % beta_hat(k) = (1 / norm( X(:,k), 'fro')^2) * soft_thresholding( r' * X(:,k) , n_dimensions*lambda );
            % MULTISAMPLE :: r' * X(:,k) is size [n_samples * 1] <-- [n_samples *  n_dimensions] * [n_dimensions*1]
            %                Should we set beta_hat(k) = mean() this result?
            %                Or should we do soft_thresholding(  mean(r'*X(:,k)) , n_dimensions*lambda)
            beta_k_tmp = (1 / norm( X(:,k), 'fro')^2) ...
                         * soft_thresholding( r' * X(:,k) , n_dimensions*lambda );
            beta_hat(k) = mean(beta_k_tmp); % AVERAGE?
            
        end % k
        
        % intercept_hat <- sum(y-X%*%beta_hat)/(n_dimensions);
        % ONE__SAMPLE
        %intercept_hat = sum(Y-X*beta_hat)) / n_dimensions;
        % MULTISAMPLE
        error('FIX intercept_hat HERE')
        intercept_hat = sum(mean(bsxfun(@minus,Y,X*beta_hat))) / n_dimensions;
        
        %intercept_list[j+1] <- intercept_hat
        intercept_hat_list(j+1) = intercept_hat;
        
        % betalist[[j+1]] <- beta_hat
        beta_list{j+1} = beta_hat;
        
        % obj[j] <- (1/2)*(1/n_dimensions)
        %         * norm(y - X%*%beta_hat - intercept_hat*rep(1,n_dimensions),"F")^2
        %         + lambda*sum(abs(beta_hat))
        % ONE__SAMPLE ::
        % obj(j) = 1/2 *( 1/n_dimensions)                                 ...
        %    * norm( y - X*beta_hat - intercept_hat*ones(n_dimensions,1),'fro')^2 ...
        %    + lambda * sum(abs(beta_hat));
        % MULTISAMPLE :: SCALE norm(Y-X*beta_hat,'fro') by sqrt(n_samples)
        obj(j) = 1/2 *( 1/n_dimensions)                                 ...
            * (norm(bsxfun(@minus,Y,X*beta_hat - intercept_hat*ones(n_dimensions,1)),'fro')/sqrt(n_samples))^2 ... % I guess I'm assuming that frob. norm([a b]) is ~ mean(norm(a),norm(b))
            + lambda * sum(abs(beta_hat));
        
        
        % if (norm(rbind(intercept_list[j],betalist[[j]]) - rbind(intercept_hat,beta_hat),"F") < tol) { break }
        % rbind should rectcle intercept_list(j) into a column vector with the column vector betalist{j}
        % Todo :: check this in R to make sure that rbind is doing what you think.
        tmp_a = [ intercept_hat_list(j) * ones(n_predictors,1)  beta_list{j} ];
        tmp_b = [ intercept_hat * ones(n_predictors,1) beta_hat ];
        if norm( tmp_a - tmp_b, 'fro') < tol && allow_break
            if verbose,fprintf('Breaking at j = %.0f\n',j),end
            break
            %             if lasso_kkt_check(X,Y,beta_hat,lambda,use_intercept,intercept_hat,verbose);
            %                 break
            %             else
            %                 if verbose(disp('   Continuing after failed check.')),end
            %                 continue
            %             end
        end
        
    end %j
    
    %warning('lasso_kkt_check not implemented')
    check = lasso_kkt_check(X,Y,beta_hat,lambda,use_intercept,[],[],intercept_hat);
    
else 
    %% use_intercept = false;
    for j = 1:maxiter
        
        for k = 1:n_predictors
            
            %r <- y - X[,-k]%*%beta_hat[-k]
            not_k = [1:(k-1) (k+1):n_predictors];
            % ONE__SAMPLE ::
            % r = y - X(:,not_k) * beta_hat(not_k);
            % MULTISAMPLE
            r = bsxfun(@minus,Y,X(:,not_k) * beta_hat(not_k));
            
            % beta_hat[k] <- (1/norm(as.matrix(X[,k]),"F")^2)*soft_thresholding(t(r)%*%X[,k],n_dimensions*lambda)
            % ONE__SAMPLE ::
            % beta_hat(k) = (1 / norm(X(:,k),'fro')^2)   ...
            %    * soft_thresholding( r' * X(:,k), n_dimensions*lambda);
            % MULTISAMPLE :: Average across samples?
            beta_k_tmp = (1 / norm(X(:,k),'fro')^2) * soft_thresholding( r' * X(:,k), n_dimensions*lambda);
            beta_hat(k) = mean(beta_k_tmp);                        
            assert(not(isnan(beta_hat(k))),'beta_hat(k) is a Nan!')
            
        end
        
        %   betalist[[(j+1)]] <- beta
        beta_list{j+1} = beta_hat;
        
        %   obj[j] <- (1/2)*(1/n_dimensions)
        %         * norm(y - X%*%beta_hat,"F")^2
        %         + lambda*sum(abs(beta_hat))
        % ONE__SAMPLE ::
        % obj(j) = 1/2 *(1/n_dimensions)             ...
        %     * norm( y - X*beta_hat,'fro')^2     ...
        %     + lambda * sum(abs(beta_hat));
        % MULTISAMPLE ::
        obj(j) = 1/2 * (1/n_dimensions)                                 ...
            * (norm( bsxfun(@minus,Y,X*beta_hat),'fro')/sqrt(n_samples))^2  ...
            + lambda * sum(abs(beta_hat));
        
        % If the ß isnt changing, break.
        % if (norm(betalist[[j]] - beta,"F") < tol) { break }
        if norm( beta_list{j} - beta_hat, 'fro') < tol && allow_break
            if verbose,
                fprintf('Breaking at j = %.0f\n',j)
            end
            break
            %             if  lasso_kkt_check(X,Y,beta_hat,lambda,use_intercept,[],[],verbose)
            %                 break
            %             else
            %                 if verbose, disp('   Continuing after failed check.'), end
            %                 continue
            %             end
        end
    end %j
    
    %warning('lasso_kkt_check not implemented')
    check = lasso_kkt_check(X,Y,beta_hat,lambda,use_intercept,[],[],verbose);

end % not(intercept_hat)


% Output.
if verbose
    if check == 1,
        fprintf('Minimum obtained.\n')
    else
        fprintf('*** Minimum not obtained. ***\n')
    end
end

beta_struct = struct();
beta_struct.beta_hat = beta_hat;

if use_intercept,
    beta_struct.intercept_hat = intercept_hat;
end

if verbose  
    fprintf('\n::: RESULTS :::\n')
    fprintf('LAMBDA :: % 2.2f\n',lambda)
    if generate_synthetic_dataset,
        fprintf('BETA\n')
        fprintf('\tTRUE\t::   RECOVERED\t::\tERROR\t::\tERROR%%\n')
        for i = 1:length(beta_true)
            fprintf('\t% 1.2f\t::\t% 1.2f\t::\t% 1.2f\t::\t% 1.1f%% \n',beta_true(i),beta_hat(i),beta_true(i)-beta_hat(i), 100*(beta_true(i)-beta_hat(i))/range(beta_true))
        end
        fprintf('\t(Error as %% of range of true beta_hat coefficients)\n')
        if use_intercept
            fprintf('intercept_hat\n')
            fprintf('\tTRUE\t::\tRECOVERED\n')
            fprintf('\t% 1.2f\t::\t% 1.2f\n',intercept,intercept_hat)
        end
    end
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = soft_thresholding(x,a)
result = zeros(size(x));
result(x >  a) = x(x >  a) - a;
result(x < -a) = x(x < -a) + a;
end


function pass = lasso_kkt_check(X,Y,beta_hat,lambda,use_intercept,intercept_hat,tol,verbose)
% See tests described at
% http://jocelynchi.com/a-coordinate-descent-algorithm-for-the-lasso-problem/#an-r-implementation-of-the-coordinate-descent-algorithm-for-the-lasso-problem
% Original code in R
%
% # Function for checking KKT conditions
% lasso_kkt_check <- function(X,y,beta_hat,lambda,intercept_hat=FALSE,intercept_hat=NULL,tol=1e-3){
%   beta_hat <- as.matrix(beta_hat); X <- as.matrix(X)
%   if (intercept_hat==TRUE){
%     Xl <- cbind(rep(1, n_dimensions),X)
%     betal <- rbind(intercept_hat,beta_hat)
%     G <- t(Xl)%*%(y-Xl%*%betal); G <- G/n_dimensions
%     if (intercept_hat - sum(y-X%*%beta_hat)/n_dimensions > tol) { return(pass=0) }
%     ix <- which(betal[2:length(betal)] == 0 )
%     iy <- which(betal[2:length(betal)] != 0)
%     if (any(abs(G[ix+1]) - lambda > tol) ) { return(pass=0) }
%     if (any(abs(G[iy+1] - lambda*sign(betal[iy+1])) > tol)) { return(pass=0) }
%   }
%   else{
%     G <- t(X)%*%(y-X%*%beta_hat)/n_dimensions
%     ix <- which(beta_hat == 0 )
%     iy <- which(beta_hat != 0)
%     if (any(abs(G[ix]) > (lambda + tol) )) { return(pass=0) }
%     if (any(G[iy] - lambda*sign(beta_hat[iy]) > tol)) { return(pass=0) }
%   }
%   return(pass=1)
% }
%

if nargin < 7 || isempty(tol)
    tol = 1E-3;
end

beta_hat = beta_hat(:);
assert(size(X,2) == length(beta_hat),'Predictor # mismatch :: # Columns of X must equal length(beta_hat)')
n_dimensions = size(X,1);
n_samples = size(Y,2);
if use_intercept
    
    Xl    = [ ones(n_dimensions,1) X];
    betal = [ intercept_hat ; beta_hat];
    % ONE__SAMPLE :: G is size [n_predictors+1 1] <----- Xl'[n_predictors+1 n_dimensions] * ( y[N_dimensions 1] - Xl[N_dimensions n_predictors+1] * betal[n_predictors+1 1])
    %G = Xl' * (y - Xl * betal);
    % MULTISAMPLE ::
    % G is size [n_predictors+1 n_samples] <------ Xl'[n_predictors+1 n_dimensions] * ( y[N_dimensions n_samples] etc.)
    G = Xl' * bsxfun(@minus, Y , Xl*betal );
    G = mean(G,2); % Average across samples ????
    G = G / n_dimensions;
    
    %  if (intercept_hat - sum(y-X%*%beta_hat)/n_dimensions > tol) { return(pass=0) }
    % ONE__SAMPLE
    % if  ((intercept_hat - sum(y - X*beta_hat)/n_dimensions) > tol)
    % MULTISAMPLE
    keyboard
    if  (mean(intercept_hat - sum( bsxfun(@minus,Y,X*beta_hat))/(n_dimensions)) > tol)
        if verbose, fprintf('lasso_kkt_check failed at intercept_hat(1) :: '); end
        pass = 0;
        return
    end
    
    ix = find( betal(2:end) == 0 ); % Coefficients == 0
    iy = find( betal(2:end) ~= 0 ); % Coefficients ~= 0
    
    %if (any(abs(G[ix+1]) - lambda > tol) ) { return(pass=0) }
    if any(  (abs(G(ix+1)) - lambda) > tol)
        if verbose, fprintf('lasso_kkt_check failed at intercept_hat(2) :: '), end
        pass = 0;
        return
    end
    
    % if (any(  abs(G[iy+1] - lambda*sign(betal[iy+1])) > tol)) { return(pass=0) }
    if   any(   abs(G(iy+1) - lambda*sign(betal(iy+1))) > tol)
        if verbose, fprintf('lasso_kkt_check failed at intercept_hat(3) :: '), end
        pass = 0;
        return
    end

else % no intercept
   
    % R :: G <- t(X)%*%(y-X%*%beta_hat)/length(y)
    % ONE__VARIABLE :: G = X' * (y - X*beta_hat)/length(y)
    % MUTLI_VARIABLE
    % Size : [n_predictors n_samples] <-- [n_predictors n_dimension] * [n_dimensions n_samples]
    % G(i,j) is the gradient(?) for (predictor_i * sample_j).
    G = X' * bsxfun(@minus,Y,X*beta_hat)/n_dimensions;
    [ind_beta_zero]   = find(beta_hat == 0); % Coefficients that have collapsed to zero
    [ind_beta_nonzero]= find(beta_hat ~= 0); % Coefficients not at 0.

    % Stationarity test.  
    % From  grad(lagrange(ß)) + lambda - lambda+ = 0
    % ..    grad(lagrange(ß)) + lambda = lambda+ >= 0
    % ..    grad(lagrange(ß)) >= -lambda
    %
    % From -grad(lagrange(ß)) + lambda - lambda- = 0
    % ..   -grad(lagrange(ß)) + lambda = lambda- >= 0
    % ..   -grad(lagrange(ß)) >= -lambda
    % ..    grad(lagrange(ß)) <=  lambda
    %
    % Therefore
    %   |grad(lagrange(ß))|  <= lambda
    %
    % Thus, for coefficients that have collapsed to zero, make sure that the
    % gradient is less than lambda + tol.
    % R :: if (any(abs(G[ix]) > (lambda + tol) )) { return(pass=0) }
    % ONE__SAMPLE ::     if any(  (abs(G(beta_zero)) - lambda) >  tol) %#ok<FNDSB>
    % MULTISAMPLE :: Strict checking across all samples.
    ktt_test1 = (abs(G(ind_beta_zero,:)) - lambda) >  tol;
    if any( ktt_test1 )
        if verbose, fprintf('lasso_kkt_check failed at no_intercept condition (1) for %.0f/%.0f samples :: ',sum(ktt_test1),n_samples), end
        pass = 0;
        return
    end
    
    % Complementary slackness test.
    % ß+ > 0, and lambda > 0 --> lambda+ === 0;
    % Since  grad(lagrange(ß)) + lambda = lambda+ = 0
    %        grad(lagrange(beta_hat)) = -lambda < 0, since lambda > 0
    % Since -grad(lagrange(ß)) + lambda = lambda- >=0
    % ..     2*lambda = lambda- > 0, since lamba > 0
    % By complementary slackness, ß- = 0.
    % Hence when ß+ > 0 -->  ß- = 0 , and grad(lagrange(ß)) = -lambda.
    % (Ditto for ß-.)
    %        
    % Thus, coefficients that have not collapsed to zero, we ensure that
    % gradient ~= lambda+ for ß+ (or lambda- for ß-).
    %
    % Note : It feels a little like there's a sign error here, but I'm trusting
    % the R code over my own intuition for now.
    % 
    % R :: if (any(G[iy] - lambda*sign(beta_hat[iy]) > tol)) { return(pass=0) }
    % ONE__SAMPLE :: if any(  ( G(ind_beta_nonzero)-lambda*sign(beta_hat(ind_beta_zero)) ) > tol )
    % MULTISAMPLE :: 
    
%     n_lambdas    = 20;
%     n_predictors = 10;
%     n_dimensions = 21;
%     n_samples    = 50;

    % G is size( n_predictors, n_samples)

    ktt_test2 = bsxfun(@minus,G(ind_beta_nonzero,:), lambda*sign(beta_hat(ind_beta_nonzero))) > tol;
    if any( ktt_test2 )
        if verbose, 
            fprintf('lasso_kkt_check failed at no_intercept(2) for an average of %.2f predictors per sample (/%.0f samples) :: ',mean(sum(ktt_test2)),n_samples), 
        end
        pass = 0;
        return
    end
    
end
% Success
pass = 1;

end





