function thresholdsize = fMRIMonteCluster(volumesize,localp,oneortwotails,globalp,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity,savename)
    % thresholdsize = fMRIMonteCluster(volumesize,localp,oneortwotails,globalp,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity,savename)
    %
    % Monte Carlo simulation to determine cluster size or mass threshold to control for multiple comparisons across % voxels in a statistical parametric map
    % and across multiple statistical maps.
    %
    % thresholdsize = minimum threshold size to achieve p < globalp corrected for multiple comparisons after local voxel-wise thresholding at p < localp
    %       threshold size is number of voxels (not mm^3) if volumeormass == 'numvoxels' or else statistical mass if volumeormass == 'mass'
    %       If cannot find thresholdsize to satisfy p < globalp corrected, thresholdsize = inf. In this case, try increasing numiter to get better sampling
    %       of larger, more rare cluster sizes.
    %
    % volumesize    = [xlen,ylen,zlen] - number of voxels along x, y, and z axes.
    % localp        = p value for local (voxel-wise) thresholding
    % oneortwotails = 2 for two-tailed or 1 for one-tailed only
    % globalp       = desired global p value(s) after cluster-mass thresholding. 
    %                 Can be a vectors - eg: [0.05 0.01 0.001] - in which case returned thresholdsize is a corresponding vector.
    % fwhm          = fwhm for Gaussian smoothing (in voxels, NOT mm).
    %                 Eg: If voxel size is 4x4x4 mm and you smoothed with Gaussian FWHM = 8 mm, pass in fwhm = 2.
    %                 Pass in fwhm = 0 to not smooth.
    %                 Uses MATLAB's smooth3 function.
    %                 Currently, only isotropic smoothing implemented.
    % mask          = [] (default) to not mask or binary volume denoting "brain" (ones) vs "not brain" (zeros)
    % volumeormass  = 'numvoxels' (default) to threshold to threshold on number of voxels in cluster
    %                 'mass' to threshold on cluster mass. Cluster mass = sum of |statistics| (eg: T-statistic absolute values)
    %                  across all voxels in the cluster, after local thresholding (eg: at p < 0.01).
    % numiter       = number iterations (default 1000)
    % distribution  = 'normal' (default) standard normally-distributed random volumes
    %                 'uniform' unit uniform-distributed random volumes
    %                 't' t-distributed random volumes, df set by localdf (below)
    %                 'f' f-distributed random volumes, df's set by localdf (below), oneortwotails must be set to 1 for f-distributions
    %                 Simulating t- and f-distributions is currently much slower than normal distributions.
    %                 NOTE: Distribution choice makes no difference with no smoothing (FHWM = 0 or []). Distribution can make increasingly
    %                 large difference with larger smoothing FWHM, especially f-distribution vs. normal.
    % localdf       = ignored for distribution == 'normal' or 'uniform'
    %               = df for distribution == 't'
    %               = [df_numerator df_denom] for distribution == 'f'
    % nummaps       = number of maps to correct across, default = 1
    %                 eg: If you compute 10 different statistical maps, you can compute a threshold
    %                 that is large enough to correct for comparisons across all voxels across all maps
    %                 by setting nummaps = 10.
    % connectivity  = 6, 18, or 26 (default), determines face vs in-plane corner vs all corner neighbours
    % savename      = filename to save results and some intermediate variables used in the simulation, default '' does not save a file
    %
    % eg: For 60 x 70 x 45 volume with voxels size 3x3x3 mm and smoothing FWHM 8mm, using local thresholding of p<0.01 two-tailed, where you want to know the
    % thresholdsize to achieve corrected p<0.05, p<0.01, and p<0.001, simulating normally-distributed random volumes, you would use this command:
    % thresholdsize = fMRIMonteCluster([60,70,45],0.01,2,[0.05, 0.01, 0.001],8/3)
    %
    %
    % Matthew Robert Graham Brown, PhD
    % October 2013
    % mrbrown23 at gmail dot com
    %
    % Copyright 2013 Matthew Robert Graham Brown
    % This software is distributed under the terms of the GNU Lesser General Public License version 3.
    % 
    % This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
    % 
    % You should have received a copy of the GNU Lesser General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


    error(nargchk(5,inf,nargin))
    if nargin < 6 | isempty(mask)          , mask          = []               ; end
    if nargin < 7 | isempty(volumeormass)  , volumeormass  = 'numvoxels'      ; end
    if nargin < 8 | isempty(numiter)       , numiter       = 1000             ; end
    if nargin < 9 | isempty(distribution)  , distribution  = 'normal'         ; end
    if nargin < 10 | isempty(localdf)      , localdf       = []               ; end
    if nargin < 11 | isempty(nummaps)      , nummaps       = 1                ; end
    if nargin < 12 | isempty(connectivity) , connectivity  = 26               ; end
    if nargin < 13                         , savename      = ''               ; end

    CheckArgsValid(volumesize,localp,oneortwotails,globalp,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity,savename)
    [csallsize, csallnum, maxclustersizes] = GetSimulatedClusters(volumesize,localp,oneortwotails,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity);
    [cumdf, unique_maxcs, nmaxcs, cummaxcs, thresholdsize] = ProcessClusterResults(csallsize, csallnum, maxclustersizes, globalp);
    if ~isempty(savename) % Optional save data
        save(savename)
    end
    DisplayResults(csallsize, csallnum, cumdf, unique_maxcs, nmaxcs, cummaxcs, globalp, thresholdsize)

% --------------------
function CheckArgsValid(volumesize,localp,oneortwotails,globalp,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity,savename)
    if length(fwhm) ~= 1
        error('You must supply fwhm == a scalar value, namely FWHM in voxels (not mm) for isotropic Gaussian smoothing kernel. Anisotropic smoothing not implemented currently.')
    end
    if ~ismember(oneortwotails,[1 2])
        error('Invalid oneortwotails value')
    end
    if ~isempty(mask) & ~isequal(size(mask),volumesize)
        error('size(mask) does not match volumesize')
    end
    if ~ismember(volumeormass,{'numvoxels','mass'})
        error('volumeormass must be ''numvoxels'' or ''mass'' - invalid volumeormass: %s',volumeormass)
    end
    if ~ismember(distribution,{'uniform','normal','t','f'})
        error('Invalid distribution: %s',distribution)
    end
    if isequal(distribution,'t') & length(localdf) ~= 1
        error('With distribution == ''t'', you must supply localdf = degrees of freedom (df)')
    end
    if isequal(distribution,'f') & length(localdf) ~= 2
        error('With distribution == ''f'', you must supply localdf = [df_numerator df_denom]')
    end
    if isequal(distribution,'f') & oneortwotails == 2
        error('With distribution == ''f'', can only use oneortwotails == 1.')
    end
    if ~ismember(connectivity,[6,18,26])
        error('Invalid connectivity: %d', connectivity)
    end

% --------------------
function [csallsize, csallnum, maxclustersizes] = GetSimulatedClusters(volumesize,localp,oneortwotails,fwhm,mask,volumeormass,numiter,distribution,localdf,nummaps,connectivity)
    % Simulates volumes, finds clusters in them, returns cluster stats
    mask            = logical(mask);
    lincutoff       = 5000;
    cslin           = zeros(lincutoff,1); % cluster sizes with size <= lincutoff, linearly indexed
    csarbsize       = []; % cluster sizes with size > lincutoff
    csarbnum        = []; % numbers of cluster sizes with size > lincutoff, matches up with csarbsize
    maxclustersizes = zeros(numiter,1); % max cluster size for a given iteration. If nummaps > 1, this is max across all maps for that iteration.
    for ii = 1:numiter
        if ii == 1 | mod(ii,20)==0, fprintf('%d / %d\n',ii,numiter), end
        [binaryvol, rawvol] = GenerateVolume(volumesize,localp,oneortwotails,fwhm,mask,distribution,localdf,nummaps); % generated simulated volume(s)
        % Counting Clusters:
        maxcs = 0;
        for vi = 1:length(binaryvol)
            cc = bwconncomp(binaryvol{vi},connectivity);
            cnum = cc.NumObjects;
            for jj = 1:cnum
                if isequal(volumeormass,'numvoxels')
                    cs = length(cc.PixelIdxList{jj});
                elseif isequal(volumeormass,'mass')
                    cs = round(sum(abs(rawvol{vi}(cc.PixelIdxList{jj}))));
                end
                maxcs = max(maxcs,cs);
                if cs <= lincutoff
                    cslin(cs) = cslin(cs) + 1;
                else
                    ind = find(csarbsize==cs);
                    if ~isempty(ind)
                        csarbnum(ind) = csarbnum(ind) + 1;
                    else
                        csarbsize = [csarbsize ; cs];
                        csarbnum = [csarbnum ; 1];
                    end
                end
            end
        end
        maxclustersizes(ii) = maxcs;
    end
    % Collating across data from iterations
    ind = find(cslin>0);
    if isempty(ind)
        cslin = [0];
    else
        cslin = cslin(1:ind(end));
    end
    if ~isempty(csarbsize)
        [csarbsize,index] = sort(csarbsize);
        csarbnum          = csarbnum(index);
        csallsize         = cat(1,[1:length(cslin)]',csarbsize);
        csallnum          = cat(1,cslin,csarbnum);
    else
        csallsize = [1:length(cslin)]';
        csallnum  = cslin;
    end

% --------------------
function [cumdf, unique_maxcs, nmaxcs, cummaxcs, thresholdsize] = ProcessClusterResults(csallsize, csallnum, maxclustersizes, globalp)
    cumdf = cumsum(csallnum);
    if cumdf(end) > 0
        cumdf = cumdf / cumdf(end); % normalizing
    end

    unique_maxcs = unique(maxclustersizes);
    nmaxcs = zeros(length(unique_maxcs),1);
    for ui = 1:length(unique_maxcs)
        nmaxcs(ui) = sum(maxclustersizes==unique_maxcs(ui));
    end
    cummaxcs = cumsum(nmaxcs);
    cummaxcs = cummaxcs / cummaxcs(end);

    thresholdsize = zeros(length(globalp),1);
    for kk = 1:length(globalp)
        ind = find(cummaxcs > (1-globalp(kk)));
        if ~isempty(ind) & ind(1) < length(cummaxcs)
            thresholdsize(kk) = unique_maxcs(ind(1));
        else
            thresholdsize(kk) = inf;
        end
    end

% --------------------
function [binaryvol, rawvol] = GenerateVolume(volumesize,localp,oneortwotails,fwhm,mask,distribution,localdf,nummaps)
    % Generates random volume(s)
    binaryvol = {};
    rawvol    = {};
    for mapindex = 1:nummaps
        switch distribution
        case 'uniform'
            [tmpbinaryvol,tmprawvol] = GenerateVolume_uniform(volumesize,localp,oneortwotails,fwhm,mask);
        case 'normal'
            [tmpbinaryvol,tmprawvol] = GenerateVolume_normal(volumesize,localp,oneortwotails,fwhm,mask);
        case 't'
            [tmpbinaryvol,tmprawvol] = GenerateVolume_t(volumesize,localp,oneortwotails,localdf(1),fwhm,mask);
        case 'f'
            [tmpbinaryvol,tmprawvol] = GenerateVolume_f(volumesize,localp,localdf(1),localdf(2),fwhm,mask);
        otherwise
            error('Unrecognized distribution: %s',distribution)
        end
        binaryvol = [ binaryvol ; tmpbinaryvol(:) ];
        rawvol    = [ rawvol    ; tmprawvol(:) ];
    end

% --------------------
function [binaryvol,rawvol] = GenerateVolume_uniform(volumesize,localp,oneortwotails,fwhm,mask)
    % Generates volume sampled from uniform distribution on the unit interval.
    vol = reshape(rand(prod(volumesize),1),volumesize);
    if fwhm > 0
        vol    = GaussianSmooth(vol,fwhm); % Gaussian smoothing of uniform volume appears to make it Gaussian
        cdfvol = normcdf(vol);
    else
        cdfvol = vol;
    end
    if ~isempty(mask)
        vol(~mask) = 0;
    end
    [binaryvol,rawvol] = ThresholdVol(cdfvol,vol,oneortwotails,localp);

% --------------------
function [binaryvol,rawvol] = GenerateVolume_normal(volumesize,localp,oneortwotails,fwhm,mask)
    % Generates standard normally-distributed volume. Similar to AlphaSim.
    vol = randn(volumesize);
    % Smoothing vol
    if fwhm > 0
        vol = GaussianSmooth(vol,fwhm);
    end
    if ~isempty(mask)
        vol(~mask) = 0;
    end
    cdfvol = normcdf(vol);
    [binaryvol,rawvol] = ThresholdVol(cdfvol,vol,oneortwotails,localp);

% --------------------
function [binaryvol,rawvol] = GenerateVolume_t(volumesize,localp,oneortwotails,df,fwhm,mask)
    if fwhm > 0 % with smoothing
        vol = zeros([volumesize df]);
        for di = 1:df
            vol(:,:,:,di) = GaussianSmooth( randn(volumesize), fwhm );
        end
    else % no smoothing
        vol = randn([volumesize df]); 
    end
    vol = mean(vol,4) ./ sqrt(var(vol,[],4) / df);
    if ~isempty(mask)
        vol(~mask) = 0;
    end
    cdfvol = tcdf(vol,df);
    [binaryvol,rawvol] = ThresholdVol(cdfvol,vol,oneortwotails,localp);

% --------------------
function [binaryvol,rawvol] = GenerateVolume_f(volumesize,localp,df_numerator,df_denominator,fwhm,mask)
    if fwhm > 0 % with smoothing
        z_numerator = GaussianSmooth( randn(volumesize), fwhm).^2;
        for ii = 1:df_numerator-1
            z_numerator = z_numerator + GaussianSmooth( randn(volumesize), fwhm).^2;
        end
        z_denominator = GaussianSmooth( randn(volumesize), fwhm).^2;
        for ii = 1:df_denominator-1
            z_denominator = z_denominator + GaussianSmooth( randn(volumesize), fwhm).^2;
        end
    else % no smoothing
        z_numerator = randn(volumesize).^2;
        for ii = 1:df_numerator-1
            z_numerator = z_numerator + randn(volumesize).^2;
        end
        z_denominator = randn(volumesize).^2;
        for ii = 1:df_denominator-1
            z_denominator = z_denominator + randn(volumesize).^2;
        end
    end
    z_numerator   = z_numerator   / df_numerator; % Chi-squared distribution scaled by df
    z_denominator = z_denominator / df_denominator;
    vol           = z_numerator  ./ z_denominator; % f-distributed volume from ratio of scaled Chi-squared's
    if ~isempty(mask)
        vol(~mask) = 0;
    end
    cdfvol    = fcdf(vol,df_numerator,df_denominator);
    binaryvol = {cdfvol > (1 - localp)}; % positives only
    rawvol    = {vol};

% --------------------
function vol = GaussianSmooth(vol,fwhm)
    gaussstd = fwhm / (2*sqrt(2 * log(2))); % see Wikipedia: fwhm
    sz = ceil(fwhm*2);
    if mod(sz,2) == 0
        sz = sz + 1;
    end
    vol = RenormalizeVolume( smooth3(vol ,'gaussian', sz, gaussstd) );

% --------------------
function normvol = RenormalizeVolume(vol)
    meanvol = mean(vol(:));
    stdvol  = std(vol(:));
    normvol = (vol-meanvol)/stdvol;

% --------------------
function [binaryvol,rawvol] = ThresholdVol(cdfvol,vol,oneortwotails,localp)
    % Thresholds cdfvol based on local p-value
    if oneortwotails == 1
        binaryvol = {cdfvol > (1 - localp)}; % positives only
        rawvol    = {vol};
    elseif oneortwotails == 2
        volpos    = cdfvol > (1 - localp/2); % positives
        volneg    = cdfvol < localp/2; % negatives
        binaryvol = {volpos,volneg};
        rawvol    = {vol ; vol};
    end

% --------------------
function DisplayResults(csallsize, csallnum, cumdf, unique_maxcs, nmaxcs, cummaxcs, globalp, thresholdsize)
    fprintf('All clusters\n')
    fprintf('   Size/Mass        Num   Cum. Prop\n')
    for kk = 1:length(csallsize)
        fprintf('[      %5d  %9d       %.3f  ]\n',csallsize(kk),csallnum(kk),cumdf(kk))
    end
    fprintf('\n')
    fprintf('Maximum cluster for each iteration (used to determine thresholdsize) \n')
    fprintf('   Size/Mass        Num   Cum. Prop\n')
    for kk = 1:length(unique_maxcs)
        fprintf('[      %5d  %9d       %.3f  ]\n',unique_maxcs(kk),nmaxcs(kk),cummaxcs(kk))
    end
    fprintf('\n  Global P          Threshold Size/Mass\n')
    for kk = 1:length(thresholdsize)
        fprintf('  %f          %f\n',globalp(kk),thresholdsize(kk))
    end

