function [sol, output, z] = LOBD(sols, R, varargin)
    %% Parse inputs 
    defaultNonneg = false;     % false
    defaultOrthogonal = true; % true
    defaultUseminf = false;    % false
    defaultCGiters = 500;
    defaultSameT = false;
    defaultMaxiters = 100; 
    defaultTvReg = false;
    defaultExpt = false;
    defaultTvWeight = 1;
    defaultProjectCoeffs = true;
    defaultShowEvery = 10;
    defaultTol = 1e-20;
    
    [U, ~, ~] = svd([sols{1}, sols{2}]);
    X = U(:, 1:R);
    [U,~] = hqrd(X);
    z = generate_z(U,R);
    D1 = (X'*sols{1}(:,1))';
    D2 = (X'*sols{2}(:,1))';
    if length(sols) == 3
        D3 = (X'*sols{3}(:,1))';
    end
    T = sols{1}'*X ./ D1 + sols{2}'*X ./ D2;
    T = T / 2;
    
    variables = struct;
    variables.x = z; %X(:,1:R); %z;
    variables.x2 = X; %X(:,1:R); %z;
    variables.t = T(1:end, 1:R) ./ T(1, 1:R);
    variables.d1 = D1 .* T(1, 1:R);
    variables.d2 = D2 .* T(1, 1:R);
    if length(sols) == 3
        variables.d3 = D3 .* T(1, 1:R);
    end
    ts = linspace(0, 1, size(sols{1}, 2))';
    variables.tz = min(exp(median(log(abs(T)) ./ ts)), 1);
%     size(variables.t);
%     size(variables.x2);
%     size(variables.d1);
%     size(variables.d2);
    
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    addRequired(p, 'sols', @(cl) all(cellfun(@(x) isnumeric(x), cl)) && length(cl) >= 2)
    addRequired(p, 'R', validScalarPosNum)
    
    addOptional(p, 'initialvalues', variables);
    addOptional(p, 'tpoints', ts)
    addParameter(p, 'maxiters', defaultMaxiters, validScalarPosNum)
    addParameter(p, 'nonneg', defaultNonneg, @islogical)
    addParameter(p, 'orthogonal', defaultOrthogonal, @islogical)
    addParameter(p, 'useminf', defaultUseminf, @islogical)
    addParameter(p, 'cgiters', defaultCGiters, validScalarPosNum)
    addParameter(p, 'sameT', defaultSameT, @islogical)
    addParameter(p, 'tvreg', defaultTvReg, @islogical)
    addParameter(p, 'showevery', defaultShowEvery, validScalarPosNum)
    addParameter(p, 'projectcoeffs', defaultProjectCoeffs, @islogical)
    addParameter(p, 'tvweight', defaultTvWeight, @isnumeric)
    addParameter(p, 'tol', defaultTol, @isnumeric)
    addParameter(p, 'expt', defaultExpt, @islogical)
    parse(p, sols, R, varargin{:})
        
    %% Set up the LOBD problem
    sols = p.Results.sols;
    Ns = length(sols);
    Nx = size(sols{1},1); % number of x points
    sz = [Nx R]; % size of the orthogonal matrix 
    orthq = @(z,task) struct_orth(z,task,sz);
    
    % Extract the initial conditions for the data
    ics = cellfun(@(x) x(:, 1)', sols, 'UniformOutput', false);
    % Constants the same for all options
    model = struct; 

    % Set up the X factor
    if p.Results.orthogonal
        if isfield(p.Results.initialvalues, 'x')
            model.variables{1} = p.Results.initialvalues.x;
        else
            xtmp = p.Results.initialvalues.variables{1};
            if all([Nx R] == size(xtmp), 'all')
                [U, ~, ~] = svd(xtmp);
                X = U(:, 1:R);
                [U,~] = hqrd(X);
                z = generate_z(U,R);
                model.variables{1} = z;
            else
                model.variables{1} = p.Results.initialvalues.variables{1};
            end
        end
        model.factors{1}  = {1, orthq};
        figure
        plot(real(orthq(model.variables{1}, [])))
    else
        if isfield(p.Results.initialvalues, 'x2')
            model.variables{1} = p.Results.initialvalues.x2;
        elseif isfield(p.Results.initialvalues, 'x')
            model.variables{1} = p.Results.initialvalues.x;
        else
            model.variables{1} = p.Results.initialvalues.variables{1};
        end
        model.factors{1}  = {1};
    end
    
    % Set up the T factor
    if ~p.Results.sameT
        if isfield(p.Results.initialvalues, 'tz') && p.Results.expt
            size(p.Results.initialvalues.tz)
            model.variables{2} = p.Results.initialvalues.tz;
        elseif isfield(p.Results.initialvalues, 't')
            model.variables{2} = p.Results.initialvalues.t;
        else
            model.variables{2} = p.Results.initialvalues.variables{2};
        end
                
        if p.Results.nonneg
            model.factors{2}  = {2, @struct_nonneg};
        elseif p.Results.expt
            model.factors{2} = {2, @(z, task) struct_exp(z, task, p.Results.tpoints)};
        else
            model.factors{2}  = 2;
            figure
            plot(real(model.variables{2}))
        end
        Doffset = 2;
        % Define the joint factorization of the matrices X and Y.
        for i = 1:Ns
            model.factorizations{i}.data = sols{i};
            model.factorizations{i}.cpd  = {1, 2, Doffset + i};
        end
    else
        if length(fieldnames(p.Results.initialvalues)) > Ns + 2 || length(fieldnames(p.Results.initialvalues)) == Ns + 1
            for i = 1:Ns 
                if isfield(p.Results.initialvalues, ['t', num2str(i)])
                    model.variables{i + 1} = p.Results.initialvalues.(['t', num2str(i)]);
                else
                    model.variables{i + 1} = p.Results.initialvalues{i + 1};
                end
            end
        else
            for i = 1:Ns 
                if isfield(p.Results.initialvalues, 't')
                    model.variables{i + 1} = p.Results.initialvalues.t;
                elseif isfield(p.Results.initialvalues, 'tz') && p.Results.expt
                    model.variables{2} = p.Results.initialvalues.tz;
                else
                    model.variables{i + 1} = p.Results.initialvalues{2};
                end
            end
        end
        
        if p.Results.nonneg
            for i = 1:Ns
                model.factors{i + 1}  = {i + 1, @struct_nonneg};
            end
        elseif p.Results.expt
            for i = 1:Ns
                model.factors{i + 1} = {i + 1, @(z, task) struct_exp(z, task, p.Results.tpoints)};
            end
        else
            for i = 1:Ns
                model.factors{i + 1}  = {i + 1};
            end
        end
        Doffset = Ns + 1;
        % Define the joint factorization of the matrices X and Y.
        for i = 1:Ns
            model.factorizations{i}.data = sols{i};
            model.factorizations{i}.cpd  = {1, i + 1, Doffset + i};
        end
    end
    
    % Set up the diagonal terms
    if p.Results.orthogonal && p.Results.projectcoeffs
        for i = 1:Ns
            curic = ics{i};
            model.factors{i + Doffset} = {1, orthq, @(z, task) struct_matvec(z, task, curic), @(z, task) struct_conj(z, task)};
        end
    else
        if ~p.Results.orthogonal && p.Results.projectcoeffs
            error('Projecting coefficients only available is X is orthogonal')
        end
        for i = 1:Ns
            if isfield(p.Results.initialvalues, ['d', num2str(i)])
                model.variables{i + Doffset} = p.Results.initialvalues.(['d', num2str(i)]);
            else
                model.variables{i + Doffset} = p.Results.initialvalues.variables{i + Doffset};
            end
            model.factors{i + Doffset} = i + Doffset;
        end
    end
    
    % Regularization terms
    if p.Results.tvreg
        for i = 1:Doffset
            model.factors{Doffset + Ns + i} = [model.factors{i}, {@(z, task) struct_fd(z, task, 1, 1)}];
            model.factorizations{Ns + i}.regL2 = Doffset + Ns + i;
            model.factorizations{Ns + i}.relweight = p.Results.tvweight;
        end
    end
    
    % Set up the optimization
    % Options for the nls solve 
    options.Display   = p.Results.showevery;   % View convergence progress every 10 iterations.
    options.TolFun    = p.Results.tol; % Stop earlier.
    options.TolX      = p.Results.tol; % Stop earlier.
    options.CGMaxIter = p.Results.cgiters;  % Recommended if structure/coupling is imposed
    options.MaxIter   = p.Results.maxiters; % Needed if there is noise! 
    
    sdf_check(model);
    
    if p.Results.useminf == 1
        %options.Algorithm = @minf_lbfgs; %@minf_ncg;
        % Solve the SDF problem.
        options.Algorithm = @minf_lbfgsdl;
        %options.Delta = 1e-10;
        %options.M = 500;
        [sol, output] = sdf_minf(model,options);
    else
        %options.Algorithm = @nls_gncgs;
         [sol, output] = sdf_nls(model,options);
    end
end

%%% Helper functions to find the Householder reflectors
function [U,R] = hqrd(X)  
    % Householder triangularization.  [U,R] = hqrd(X);
    % Generators of Householder reflections stored in U.
    % H_k = I - U(:,k)*U(:,k)'.
    % prod(H_m ... H_1)X = [R; 0]
    % where m = min(size(X))
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    U = zeros(size(X));
    m = min(n,p);
    R = zeros(m,m);
    for k = 1:min(n,p)
        [U(k:n,k),R(k,k)] = housegen(X(k:n,k));
        v = U(k:n,k)'*X(k:n,k+1:p);
        X(k:n,k+1:p) = X(k:n,k+1:p) - U(k:n,k)*v;
        R(k,k+1:p) = X(k,k+1:p);
    end
end

function [u,nu] = housegen(x)
    % [u,nu] = housegen(x)
    % Generate Householder reflection.
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    % [u,nu] = housegen(x).
    % H = I - uu' with Hx = -+ nu e_1
    %    returns nu = norm(x).
    u = x;
    nu = norm(x);
    if nu == 0
        u(1) = sqrt(2);
        return
    end
    u = x/nu;
    if u(1) >= 0
        u(1) = u(1) + 1;
        nu = -nu;
    else
        u(1) = u(1) - 1;
    end
    u = u/sqrt(abs(u(1)));
end

function [z] = generate_z(U,R)
    Nx = size(U,1); 
    z = zeros(Nx*R-0.5*R*(R-1),1);
    zind = 1;
    for i = R:-1:1 
        z(zind:zind+Nx-i) = U(i:end,i);
        zind = zind + 1+Nx-i;
    end
end
