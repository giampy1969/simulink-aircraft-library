function [X0,U0]=air3m(name,V,H,G)

% [X0,U0]=air3m(name,V,H,G) finds aircraft trim condition
%
% For each entry in the matrices V,H,G (speed,altitude,gamma)
% the function finds the corresponding trim point condition on 
% state and input, (i.e. alpha, theta, thrust, and elevators).
% 
% If the fourth argument, gamma = theta - alpha = asin(Hdot/V),
% (also called "flight path angle"), is zero (radians), then the 
% trim will corresponds to a straight and level flight condition,
% otherwise the trim is not really an equilibrium since the
% airplane is constantly changing its altitude.
% 
% Example 1: find the straight and level flight condition for 
% the scheme airtrim.mdl corresponding to V=260 m/s and H=100 m :
% [x0,u0]=air3m('airtrim',260,100,0);
% 
% Example 2: find the two trim conditions for the points:
% 1) V=200 m/s, H=10e3 m, descending with an angle of 20 degrees
% 2) V=150 m/s, H=50 m, raising with an angle of 45 degrees
% [X0,U0]=air3m('airtrim',[200 150],[10e3 50],[-20 45]*pi/180);
% squeeze(X0),squeeze(U0), % to see the results more clearly
% 
% Example 3: find the trim conditions for the extreme points
% in the envelope 60 m/s < V < 200 m/s, 60 m < H < 16000 m
% and 0 < G < pi/6, for the scheme airtrim.mdl :
% [V,H,G]=ndgrid([60 200],[60 8000 16000],[0 pi/6]);
% [X0,U0]=air3m('airtrim',V,H,G);
% the trim point corresponding to 60 m/s, 16000 m, G=pi/6, 
% (that is V(1,3,2), H(1,3,2), G(1,3,2)) is X0(:,:,1,3,2)
% the full matrices of alpha and elevators corresponding
% to the straight and level flight condition G=0 are:
% squeeze(X0(2,:,:,:,1)) and squeeze(U0(10,:,:,:,1))
% 
% Note that the outputs matrices are ready to be used with the 
% "interpolate matrix" block from the Aerospace Blockset library.
% The function could be easily extended to give as outputs also
% the controller matrices for each trim condition, those matrices 
% too would be ready to be used with the "3D controller" block, 
% still from the GNC sublibrary of the Aerospace Blockset.

% Giampy, May 2003
% the function heavily relies on the function jj_trim from
% the great "trimmod" utility, (thanks to J. J. Buchholz).

% max number of steps and final cost function value
opt=[10e3 1e-13];

% indices of trim variables: alpha, theta, (=> idX=[2 8]')
% thrust, and elevator (=> idU=[1 7]') , and and indices 
% for the requirements on the time derivatives of x:
% vdot=0,alphadot=0,qdot=0,Hdot=V*sin(G) (=> idR=[1 2 5 12]')
idX=[2 8]';idU=[1 7]';idR=[1 2 5 12]';

% NOTE: if using stabilators for the trim, idU=[1 10]';
% remember that an airlib-like structure is assumed
% regarding the number and the order of inputs and states

% check inputs
if nargin ~= 4, 
    error('must supply 4 arguments'); 
end

if exist(name) ~= 4, 
    error(['The model ' name ' hasn''t been found']);
end

SZS=[size(V,1) size(H,1) size(G,1); ...
     size(V,2) size(H,2) size(G,2); ...
     size(V,3) size(H,3) size(G,3)];

if any(SZS<1) | any(SZS>3),
    error(['V H and G must be non empty 2D or 3D matrices']);
end

if any(SZS(:,1)~=SZS(:,2)) | any(SZS(:,2)~=SZS(:,3)), 
    error(['V H and G must have the same size']); 
end

if max(max(max(abs(G)))) > pi/2, 
    error('Gamma must be in radians and between -pi/2 and pi/2');
end

% check model
[msizes,x0,str,ts]=feval(name,0,0,0,0);
if any(msizes([1 4])~=[12;10]),
    error('model should be airlib-like');
else
    % complete initial conditions 
    % (x0 and u0 are used for the max step size computation) 
    u0=zeros(10,1);u0(1)=0.2;
    dx0=zeros(size(x0));
    y0=zeros(msizes(3),1);
end

% initialize outputs
X0=zeros(12,1,size(V,1),size(V,2),size(V,3));
U0=zeros(10,1,size(V,1),size(V,2),size(V,3));

% default linearization step sizes
dxL=1e-5*ones(size(x0));duL=1e-5*ones(size(u0));

% default maximum step sizes
dxM=1e42*ones(size(x0));duM=1e42*ones(size(u0));

% adjust maximum step size for the 4 trim variables
dxM(idX)=0.001+abs(x0(idX)/10);duM(idU)=0.001+abs(u0(idU)/10);
            
% names required by jj_3m
nx=length(x0);xnm=cell(nx,1);for i=1:nx,xnm{i}=['x_' num2str(i) ' ']; end
nu=length(u0);unm=cell(nu,1);for i=1:nu,unm{i}=['u_' num2str(i) ' ']; end
ny=length(y0);ynm=cell(ny,1);for i=1:ny,ynm{i}=['y_' num2str(i) ' ']; end
nd=length(dx0);dnm=cell(nd,1);for i=1:nd,dnm{i}=['d_' num2str(i) ' ']; end

for i = 1:size(V,1),
    for j = 1:size(V,2),
        for k = 1:size(V,3),
            
            % reset initial point
            x0=zeros(12,1);u0=zeros(10,1);u0(1)=0.2;
            
            % set new position in the flight envelope
            x0([1 12])=[V(i,j,k) H(i,j,k)];
            
            % requirement for Hdot
            dx0(12)=V(i,j,k)*sin(G(i,j,k));
            
            % finds the trim value
            [x0,u0,dx0,y0]=jj_3m(name,x0(:),u0(:),dx0(:),y0(:),idX,idU,idR,[],xnm,unm,dnm,ynm,dxM(:),duM(:),dxL(:),duL(:),opt);
            
            % write the trim value into X0 and U0
            X0(:,1,i,j,k)=x0(:);U0(:,1,i,j,k)=u0(:);
            
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jj_3m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ...
  [x_tr, u_tr, d_tr, y_tr, varargout] = jj_3m (sys, ...
  x, u, d, y, ...
  i_x, i_u, i_d, i_y, ...
  x_nam, u_nam, d_nam, y_nam, ...
  del_x_max, del_u_max, del_x_lin, del_u_lin, ...
  options, ...
  x_min, x_max, u_min, u_max)
%JJ_TRIM   Trim point determination of a nonlinear ordinary differential
%equation system
%
%   [X_TR, U_TR, D_TR, Y_TR] = JJ_TRIM (SYS, X, U, D, Y, I_X, I_U, I_D, I_Y, X_NAM, U_NAM, D_NAM, Y_NAM)
%   trims the system SYS towards an operating point defined by the elements
%   of
%   X (state vector), U (input vector), D (derivative of the state vector),
%   and Y (output vector).
%   I_X and I_U define the indices of the so-called trim variables, which
%   are those states and inputs the trim algorithm has to find
%   the appropriate values for. Specified values of the trim variables
%   are taken as initial starting guesses for the iteration.
%   I_D and I_Y are the indices of the so-called trim requirements, which
%   the trim algorithm has to satisfy. The values of the other D(i) and Y(i)
%   do not matter.
%   X_NAM, U_NAM, D_NAM, and Y_NAM are cell arrays of the form
%   X_NAM = {'state_1'; 'state_2'; ...}
%   containing the names of the states, inputs, derivatives, and outputs.
%   The names can be chosen arbitrarily. They are used only to identify
%   linear dependent trim variables or trim requirements.
%
%   IMPORTANT: o There have to be as many trim variables as there are
%                trim requirements.
%
%              o All vectors (and cell arrays) have to be column vectors.
%
%   [X_TR, U_TR, D_TR, Y_TR, ERRFLG] = JJ_TRIM (SYS, ...)
%   ERRFLG = 0 if trimming was successful,
%   ERRFLG = 1 if an error occurred during trimming.
%              The returned trim result is invalid.
%
%   [X_TR, ..., ERRFLG, INFOSTRUCT] = JJ_TRIM (SYS, ...)
%   additionally returns an INFOSTRUCT structure containing
%   - x: the progression of all states,
%   - u: the progression of all inputs,
%   - d: the progression of all state derivatives,
%   - y: the progression of all outputs
%   over all iterations. They are of size N x (N_iter+1) where N is the
%   respective vector length (number of states, inputs etc.) and N_iter is
%   the number of iterations. The first column of each quantity represents
%   the initial point (before the first iteration).
%
%   - n_iter: the total number of iterations
%
%   - cost:   the progression of the cost function (length: N_iter+1)
%
%   - jaco:   the progression of the Jacobian as a cell array of N_iter
%			  cells, each cell containing the Jacobian (size N_t_v x N_t_v)
%             of the current iteration, where N_t_v is the number of trim
%			  variables / trim requirements.
%             The rows represent trim variables (first the states x, then
%             the inputs u), the columns represent the trim requirements
%             (state derivatives d followed by outputs y). The names can be
%             extracted using the following code:
%
%             x_u_nam = [x_nam; u_nam];       	% Names of all generalized inputs
%             d_y_nam = [d_nam; y_nam];      	% Names of all generalized outputs
%             i_t_v   = [i_x; length(x)+i_u]; 	% Indices of trim variables in x_u_nam
%             i_t_r   = [i_d; length(d)+i_y];   % Indices of trim requirements in d_y_nam
%             T_V_nam = x_u_nam(i_t_v);         % Names of trim variables
%             T_R_nam = d_y_nam(i_t_r);         % Names of trim requirements
%
%             Note that the Jacobian is not updated if the cost function
%             deteriorated from the previous step (StepType 2, see below).
%             In this case, an empty  matrix is stored instead.
%
%   - StepType: enumeration for the type of each step:
%         1: regular iteration
%         2: bisection to recover from a previous iteration step which lead
%            to an increased cost function
%         3: bisection to recover from a previous iteration step which lead
%            to a rank-deficient Jacobian even though the Jacobian had full
%            rank at the/a previous step
%    The i-th entry of StepType corresponds to the step from x(i) and u(i)
%    to x(i+1) and u(i+1).
%
%   - BisecCounter: number of bisections (including the current step) since
%                   the last successful iteration
%
%   - LimitedStep: flag indicating if the step size was truncated:
%     	 0: step was not truncated
%        1: step was truncated to the user-specified step size (DEL_X_MAX
%           and/or DEL_U_MAX)
%        2: step was truncated to avoid violation of upper or lower bounds
%           (X_MIN, X_MAX, U_MIN, or U_MAX)
%        3: step was truncated because of step size limitations AND hard
%           bounds.
%
%   - ExitCode: enumeration for the (reasons of) success or failure of the
%               trimming process:
%     	 0: successful trim
%        1: trim failed:  singular Jacobian at problem setup
%        2: trim failed:  max. number of bisections reached when trying to
%                         restore a regular Jacobian
%        3: trim failed:  max. number of bisections reached when trying to
%                         reduce the cost function
%        4: trim failed:  max. number of iterations reached
%        5: trim failed:  no more step possible in the required direction
%                         without violating the specified bounds
%
% 	- OutputMessage: string or cell of strings containing a detailed message
%                    about the trim result.
%
%
%   [X_TR, ...] = JJ_TRIM (..., Y_NAM, DEL_X_MAX, DEL_U_MAX)
%   allows the additional specification of maximum (absolute) step sizes of
%   states and inputs between consecutive iterations.
%   The lengths of DEL_X_MAX and DEL_U_MAX must be equal to those of
%   X and U respectively.
%   The default values of DEL_X_MAX and DEL_U_MAX are 1e42.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_MAX, DEL_X_LIN, DEL_U_LIN)
%   allows the additional specification of the state and input perturbations
%   to be used for the calculation of the Jacobian matrix (sensitivity matrix)
%   in the linearization procedure. Their length must match the respective
%   lengths of x and u. All elements must be nonzero and are automatically
%   treated as absolute values. The default values of DEL_X_LIN and DEL_U_LIN
%   are 1e-6*(1 + abs(x)) and 1e-6*(1 + abs(u)), respectively.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_LIN, OPTIONS)
%   allows the specification of additional options via designated fields of
%   the OPTIONS structure:
%   n_iter_max:   	the maximum number of iterations (default: 42)
%
%   cost_tbg:       the cost value to be gained (default: 1e-9)
%
%   CompileFlag:    flag {0,1,2} to control if/how to compile the model (default: 1)
%       o flag = 0: model is not compiled by jj_trim (hence needs to be precompiled)
%       o flag = 1: model is compiled by jj_trim using the 'compile' command.
%                   This is the default option in order to be backward
%                   compatible with jj_trim version 1.5 and lower.
%       o flag = 2: model is compiled by jj_trim using the 'lincompile'
%                   command. Recommended e.g. if the model incorporates time
%                   delay blocks to be used with direct feedthrough during
%                   linearization.
%
%   EnableMessages: flag to disable (0) or enable (~=0) message boxes and
%                   console messages. Messages are enabled by default. Note that
%                   disabling messages will make the function completely silent.
%
%	n_bisec_max:	Maximum number of consecutive bisections (default: 10)
%
%   If a field is not set (i.e. is not part of the structure), the
%   respective default value is used.
%   In order to maintain backward compatibility with jj_trim versions 1.5
%   and below, it is also possible to define OPTIONS as a 2-element vector.
%   In this case the following definition applies:
%	- OPTIONS(1): the maximum number of iterations (default: 42)
%   - OPTIONS(2): the cost value to be gained (default: 1e-9)
%   If OPTIONS is passed as a vector, it is not possible to set CompileFlag,
%   EnableMessages, or n_bisec_max, and their default values will be used.
%
%   OPTIONS can also be passed as an empty vector [], in which case the
%   default values will be used.
%
%
%   [X_TR, ...] = JJ_TRIM (..., OPTIONS, X_MIN, X_MAX, U_MIN, U_MAX)
%   allows the additional specification of fixed upper and lower boundaries
%   for states x and inputs u. They can be defined as vectors, in this case
%   their lengths must match the respective lengths of x and u.
%   You can use an empty vector [] to use the default values: -inf for
%   X_MIN and U_MIN, and inf for X_MAX and U_MAX, respectively.
%
%
%   Function dependencies:
%   1. jj_lin (for linearization)
%   2. getGeneralizedOutputVector (model evaluation)
%
%
%   Notes:
%   Unfortunately, Simulink claims the right to alter the order of the states
%   during the simulation (see https://mathworks.com/help/simulink/gui/initial-state.html).
%   MathWorks hence recommends not to use arrays for initialization. There
%   are different solutions to circumvent this; one solution is converting
%   the array to a structure (additionally features the corresponding state
%   names) and initializing the model using this structure as shown below:
%
%   u_tr_with_zero = [0, u_tr']; % add a zero as the time span trailing the input vector
%   x_tr_string = mat2str (x_tr, 42); % convert the state trim vector to a string with maximum precision:16
% 
%   set_param(sys,'LoadInitialState','on');  % enable the checkboxes in the Initial state...
%   set_param(sys,'LoadExternalInput','on'); % and Input menu entry of the current model
%   set_param(sys,'InitialState',x_tr_string ); % transfer the state trim vector to the corresponding edit text
%   set_param(sys,'SaveFormat','Structure'); % set the SaveFormat to Structure
% 
%   x_tr_structure = Simulink.BlockDiagram.getInitialState(sys); % reread the initial state back as a structure into a new variable
%   model_workspace = get_param(sys,'ModelWorkspace'); % retrieve a handle to the model workspace
% 
%   assignin(model_workspace,'initial_state_in_model_workspace',x_tr_structure); % write the initial state (as a structure)...
%   assignin(model_workspace,'initial_input_in_model_workspace',u_tr_with_zero); % and the initial input (as a vector) to the model
% 
%   set_param(app.model_name,'ExternalInput','initial_input_in_model_workspace'); % use these model workspace variables...
%   set_param(app.model_name,'InitialState','initial_state_in_model_workspace');  % in the Input and Initial state edit texts
%   Authors:
%   J. J. Buchholz, Hochschule Bremen,              buchholz@hs-bremen.de
%   W. Moennich,  	German Aerospace Center (DLR),	wulf.moennich@dlr.de
%   D. Kiehn,       German Aerospace Center (DLR),	daniel.kiehn@dlr.de
%   Version 2.2    25.02.2022, D. Kiehn
%   - Fixed a bug which caused the model not to be terminated in case the
%     trim failed due to hard bounds (i.e., exit code 5)
%   - It is now verified that the target cost function (cost_tbg) is non-negative.
%   - Added a note about the risks of directly using the trim result (as an array)
%     for the initialization of a model
%   Version 2.1    05.08.2021, D. Kiehn
%   - Fixed a bug which caused the step direction to change sign in case
%     the step had to be truncated to avoid violating the lower bounds
%     (X_MIN or U_MIN).
%   - OPTIONS can now be passed as an empty element.
%   Version 2.0     22.07.2020, D. Kiehn
%
%   - Individual checking of max. iteration step sizes and max. perturbations
%     for linearization.
%   - The suppression of the array initialization warning is now done via
%     try/catch to ensure compatibility with older Matlab versions.
%   - Added automatic bisection if the sensitivity matrix loses full rank
%     (but was initialized with full rank) - useful in case a saturation is
%     reached during the trimming process.
%   - The OPTIONS input is now a structure which can be used to
%       - control if/how to compile the Simulink model, granting the
%         ability to use precompiled models,
%       - define the maximum number of iterations and bisections,
%       - enable/disable all messages and warnings.
%   - Added optional inputs to specify fixed upper and lower boundaries for
%     states and inputs: X_MIN, X_MAX, U_MIN, U_MAX.
%   - Added optional output INFOSTRUCT with details about the trimming
%     process (progression of state/input vectors, cost function, Jacobian,
%     etc. over all iterations).
%   - Introduced a "fail-fast" philosophy: the function now returns an
%     error if it is called with nonsensical inputs, e.g. if the number of
%     trim variables doesn't match the number of trim requirements.
%   - Moved the model evaluation to new function getGeneralizedOutputVector
%
%   Note: default options of version 2.0 are those of jj_trim version 1.5
%         to maintain backward compatibility.
%   Version 1.5     07.11.2011, J. J. Buchholz
%
%   Array initialization warning suppressed.
%   Minor code optimizations.
%   Version 1.2.4   15.10.2008, W. Moennich
%
%   Comment added:
%   % The functionalities used here ('compile', 'term', ...)
%   % are documented in the Simulink manual under Simulation Commands/model
%   % and in Help(doc) unter Simulink/Functions/Simulation/model.
%   Version 1.2.3   20.05.2003, W. Moennich
%
%   Optional error flag as 5th output.
%   Version 1.2.2   26.09.2001, W. Moennich
%
%   Empty vectors now possible for DEL_X_MAX, DEL_U_MAX, DEL_X_LIN,
%   DEL_U_LIN.
%   Version 1.2.1   29.08.2001, W. Moennich
%
%   Check for maximum number of iterations moved inside the loop.
%   Version 1.2 	26.05.2000, J. J. Buchholz
%
%   The names of inputs, ... outputs for the error messages
%   are now transferred via the parameter list.
%
%   The precompilation and the release of the system is now done in JJ_TRIM.
%   Version 1.0     1998,       J. J. Buchholz
%
%   Initial version.
%% Initialize
% Open the system block diagram internally without bringing the system
% block diagram to front
dummy = eval (sys);
% Suppress the array initialization warning (only required in newer Matlab
% versions and will crash in older Matlab versions, hence the try/catch)
try
	set_param (sys, 'InitInArrayFormatMsg', 'None');
catch
end
% Feed through all initial vectors, useful in case of an emergency exit
x_tr = x;           % States
u_tr = u;           % Inputs
d_tr = d;           % State derivatives
y_tr = y;           % Outputs
varargout{1} = [];  % Error flag
varargout{2} = {};  % Infostruct
% Determine lengths of basic vectors
n_x  = length (x);
n_u  = length (u);
n_d  = length (d);
n_y  = length (y);
n_i_x = length (i_x); % Number of trim variables (states)
n_i_u = length (i_u); % Number of trim variables (inputs)
n_i_d = length (i_d); % Number of trim requirements (state derivatives)
n_i_y = length (i_y); % Number of trim requirements (outputs)
% Assemble generalized input vector and generalized output vector
x_u = [x; u];
d_y = [d; y];
% Determine length of generalized input vector
n_x_u = n_x + n_u;
% Assemble trim variable index vector and trim requirement index vector
i_t_v = [i_x; i_u + n_x];
i_t_r = [i_d; i_y + n_d];
% Determine number of trim variables and trim requirements
n_t_v = n_i_x + n_i_u;
n_t_r = n_i_y + n_i_d;
% Set evaluation time to zero explicitly (assume a time invariant system)
t = 0;
%% Read and check inputs
% Options have to be checked first to decide whether or not to display
% error messages
% If no options have been defined,
if nargin < 18
	% Set defaults which match jj_trim 1.5 and lower
	n_iter			= 42;	% Max. number of iterations (including bisection steps)
	cost_tbg		= 1e-9;	% Cost function to be gained
	CompileFlag		= 1;	% Controls if/how to compile the Simulink model (see function help)
	EnableMessages	= 1;	% Disables (0) or enables (~=0) pop-up message boxes
	N_bisec_max     = 10;	% Max. number of consecutive bisections
else
	if isempty(options) % if options are passed as empty vector, use defaults
		n_iter = 42;
		cost_tbg = 1e-9;
		CompileFlag = 1;
		EnableMessages = 1;
		N_bisec_max = 10;
	elseif isstruct(options)
		% If the options are defined as a struct, read the specified options
		% individually
		% Read max. number of iterations
		if isfield(options,'n_iter_max')
			n_iter = options.n_iter_max;
		else
			n_iter = 42;
		end
		% Read cost value to be gained
		if isfield(options,'cost_tbg')
			if options.cost_tbg >= 0
				cost_tbg = options.cost_tbg;
			else
				errordlg ('The cost value to be gained must not be negative!','Error');
				return
			end
		else
			cost_tbg = 1e-9;
		end
		% Read the compile flag
		if isfield(options,'CompileFlag')
			CompileFlag = options.CompileFlag;
		else
			CompileFlag = 1;
		end
		% Read the flag that enables/disables messages
		if isfield(options,'EnableMessages')
			EnableMessages = options.EnableMessages;
		else
			EnableMessages = 1;
		end
		% Read the max. number of consecutive bisections
		if isfield(options,'n_bisec_max')
			N_bisec_max = options.n_bisec_max;
		else
			N_bisec_max = 10;
		end
	else
		% If the "options" input is not a struct, try to access the vector
		% entries. This is required to maintain backward compatibility
		% with versions 1.5 and below.
		n_iter   = options(1);
		cost_tbg = options(2);
		EnableMessages = 1;
		CompileFlag	   = 1;
		N_bisec_max = 10;
	end
end
% Check if there are as many trim variables as there are trim requirements
if n_t_r ~= n_t_v
	l1 = ['The number of trim variables: ', int2str(n_t_v)];
	l2 = 'does not equal';
	l3 = ['the number of trim requirements: ', int2str(n_t_r)];
	l4 = ' ';
	l5 = 'The returned trim point is not valid.';
	outputMessage = {l1; l2; l3; l4; l5};
	errordlg (outputMessage,'Error');
	return
end
% If there are no trim variables and requirements
if ~n_t_r
	if EnableMessages
		% Prepare output message
		outputMessage = {'There should be at least';
						 'one trim variable and';
						 'one trim requirement.'};
		warndlg(outputMessage,'Nothing to trim');
	end
	return
end
% If no maximum step sizes have been defined for x and u
if nargin < 14
	% set defaults
	del_max = 1.e42*ones (n_x_u, 1);
	% If the maximum step sizes have been defined at least for x
else
	% Check del_x_max, set to default if empty
	if isempty(del_x_max)
		del_x_max = 1.e42*ones (n_x,1);
	elseif numel(del_x_max) ~= n_x
		outputMessage = {...
			['The number of elements in del_x_max (',num2str(numel(del_x_max)),')'];
			['must match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end
	% If no maximum step sizes have been defined for u,
	if nargin < 15
		del_u_max = 1.e42*ones (n_u, 1);
	else
		% Check del_u_max, set to default if empty
		if isempty(del_u_max)
			del_u_max = 1.e42*ones (n_u, 1);
		elseif numel(del_u_max) ~= n_u
			outputMessage = {...
				['The number of elements in del_u_max (',num2str(numel(del_u_max)),')'];
				['must match the length of u (' num2str(n_u),').']};
			errordlg (outputMessage,'Error');
			return
		end
	end
	% Assemble generalized maximum step vector
	del_max = [del_x_max(:); del_u_max(:)];
end
% If no step sizes for the linearization have been defined,
if nargin < 16
	del_lin = 1e-6*(1 + abs(x_u)); % set defaults
else
	% Check del_x_lin, set to default if empty
	if isempty(del_x_lin)
		del_x_lin = 1e-6*(1 + abs(x));
	elseif numel(del_x_lin) ~= n_x
		outputMessage = {...
			['The number of elements in del_x_lin (',num2str(numel(del_x_lin)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg (outputMessage,'Error');
		return
	end
	% If no maximum step sizes have been defined for u,
	if nargin < 17
		del_u_lin = 1e-6*(1 + abs(u)); % set defaults
	else
		% Check del_u_lin, set to default if empty
		if isempty(del_u_lin)
			del_u_lin = 1e-6*(1 + abs(u));
		elseif numel(del_u_lin) ~= n_u
			outputMessage = {...
				['The number of elements in del_u_lin (',num2str(numel(del_u_lin)),')'];
				['does not match the length of u (' num2str(n_u),').']};
			errordlg(outputMessage,'Error');
			return
		end
	end
	% Assemble generalized perturbation vector for linearization
	% Make sure only positive values are used
	del_lin = [abs(del_x_lin(:)); abs(del_u_lin(:))];
end
% Check upper and lower bounds for x and u
if nargin < 19
	% If no bounds were specified at all, set defaults
	x_min = -inf*ones(n_x,1);
	x_max =  inf*ones(n_x,1);
	u_min = -inf*ones(n_u,1);
	u_max =  inf*ones(n_u,1);
else
	% If at least x_min was specified,
	if isempty(x_min)
		x_min = -inf*ones(n_x,1);
	elseif numel(x_min) == n_x
		x_min = x_min(:);
	else
		outputMessage = {...
			['The number of elements in x_min (',num2str(numel(x_min)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end
	% Check conformance of x_max with initial x
	if any(x(:) < x_min(:))
		errordlg('The initial state vector x violates the bounds of x_min!','Error');
		return
	end
end
% Check x_max
if nargin > 19
	if isempty(x_max)
		x_max = inf*ones(n_x,1);
	elseif numel(x_max) == n_x
		x_max = x_max(:);
	else
		outputMessage = {...
			['The number of elements in x_max (', num2str(numel(x_max)),')'];
			['does not match the length of x (' num2str(n_x),').']};
		errordlg(outputMessage,'Error');
		return
	end
	% Check conformance of x_max with initial x
	if any(x(:) > x_max(:))
		errordlg('The initial state vector x violates the bounds of x_max!','Error');
		return
	end
else
	x_max =  inf*ones(n_x,1);
end
% Check u_min
if nargin > 20
	if isempty(u_min)
		u_min = -inf*ones(n_u,1);
	elseif numel(u_min) == n_u
		u_min = u_min(:);
	else
		outputMessage = {...
			['The number of elements in u_min (',num2str(numel(u_min)),')'];
			['does not match the length of u (' num2str(n_u),').']};
		errordlg(outputMessage,'Error');
		return
	end
	% Check conformance of u_min with initial u
	if any(u < u_min)
		errordlg('The initial input vector u violates the bounds of u_min!','Error');
		return
	end
else
	u_min = -inf*ones(n_u,1);
end
% Check u_max
if nargin > 21
	if isempty(u_max)
		u_max = inf*ones(n_u,1);
	elseif numel(u_max) == n_u
		u_max = u_max(:);
	else
		outputMessage = {...
			['The number of elements in u_max (', num2str(numel(u_max)),')'];
			['does not match the length of u (' num2str(n_u),').']};
		errordlg(outputMessage,'Error');
		return
	end
	% Check conformance of u_max with initial u
	if any(u > u_max)
		errordlg('The initial input vector u violates the bounds of u_max!','Error');
		return
	end
else
	u_max = inf*ones(n_u,1);
end
x_u_min = [x_min; u_min];
x_u_max = [x_max; u_max];
%% Core algorithm
% Pre-allocate infostruct (only required if 6th output is desired)
if nargout > 5
	infostruct.x    = zeros(n_x,n_iter+1);
	infostruct.u    = zeros(n_u,n_iter+1);
	infostruct.d    = zeros(n_d,n_iter+1);
	infostruct.y    = zeros(n_y,n_iter+1);
	infostruct.cost = zeros(1,n_iter+1);
	infostruct.jaco = cell(1,n_iter);
	infostruct.StepType     = zeros(1,n_iter);
	infostruct.BisecCounter = zeros(1,n_iter);
	infostruct.LimitedStep  = zeros(1,n_iter);
	infostruct.n_iter   = 0;
	infostruct.ExitCode = 0;
end
% Set old cost value to infinity, in oder to definitely have 
% an improvement with the first try
cost_old = inf;
% Precompile system.
% Unfortunately, this is necessary because only precompiled systems can be evaluated.
% If the trim algorithm is aborted without the corresponding "Release system" command
% the next precompilation attempt will lead to an error and the simulation cannot
% be started.
% The system then has to be released manually (maybe more than once!) with:
% model_name ([], [], [], 'term')
% The functionalities used here ('compile', 'term', ...) 
% are documented in the Simulink manual under Simulation Commands/model
% and in Help(doc) unter Simulink/Functions/Simulation/model.
if CompileFlag == 1
	feval (sys, [], [], [], 'compile');
elseif CompileFlag == 2
	feval (sys, [], [], [], 'lincompile');
end
% Initialize bisection counter
i_bisec = 0;
% Loop over maximum n_iter iterations
for i_iter = 0 : n_iter
	% Calculate generalized output vector at the current trim point.
	d_y_tr = getGeneralizedOutputVector(sys,x_u,n_x,t);
	% Update function outputs
	x_tr = x_u(1:n_x,1);
	u_tr = x_u(n_x+1:end,1);
	d_tr = d_y_tr(1:n_d,1);
	y_tr = d_y_tr(n_d+1:end,1);
	% Calculate differences between required and current generalized output vectors
	del_d_y = d_y - d_y_tr;
	% Pick trim requirements only
	del_t_r = del_d_y(i_t_r);
	% Cost value is the maximum element of the trim requirement error vector
	cost = max(abs(del_t_r));
	% Update infostruct (if output is desired)
	if nargout > 5
		infostruct.x(1:n_x,i_iter+1) = x_u(1:n_x,1);
		infostruct.u(1:n_u,i_iter+1) = x_u(n_x+1:end,1);
		infostruct.d(1:n_d,i_iter+1) = d_y_tr(1:n_d,1);
		infostruct.y(1:n_y,i_iter+1) = d_y_tr(n_d+1:end,1);
		infostruct.cost(1,i_iter+1) = cost;
	end
	% If current cost value has become smaller 
	% than the cost value to be gained
	if cost < cost_tbg
		% Assemble output message
		l1 = ['A cost value of ',num2str(cost)];
		l2 = ['has been gained after ', int2str(i_iter), ' iteration(s).'];
		h1 = 'Success';
		outputMessage = [l1,' ',l2];
		% Display success message in console
		if EnableMessages
            helpdlg({l1, l2}, h1);
		end
		% Release system
		if (CompileFlag == 1) || (CompileFlag == 2)
			feval (sys, [], [], [], 'term');
		end
		% If a 5th output parameter is used,
		if nargout > 4
			% Return errorflag as 5th parameter
			errflg = 0;
			varargout{1} = errflg;
			% If a 6th output parameter is used,
			if nargout > 5
				% Update infostruct and remove unassigned parts
				infostruct.ExitCode	= 0; % Exit code 0: successful trim
				infostruct.OutputMessage = outputMessage;
				infostruct = truncateInfoStruct(infostruct,i_iter);
				% Assign infostruct to second output
				varargout{2} = infostruct;
			end
		end
		% Trimming was successful
		return
	end
	% If the maximum number of iterations is reached,
	% output error message and abort program
	if i_iter == n_iter
		l1 = ['Maximum number of iterations reached: ', int2str(i_iter)];
		l2 = 'Program was aborted';
		l3 = ['with a cost value of: ', num2str(cost)];
		outputMessage = {l1; l2; l3};
		if EnableMessages
			errordlg (outputMessage, 'Program aborted');
		end
		% Release system
		if (CompileFlag == 1) || (CompileFlag == 2)
			feval (sys, [], [], [], 'term');
		end
		% If a 5th output parameter is used,
		if nargout > 4
			% Update error flag and return it as 5th parameter
			errflg = 1;
			varargout{1} = errflg;
			% If a 6th output parameter is used,
			if nargout > 5
				% Update infostruct and remove unnecessary fields
				infostruct.ExitCode	= 4; % Exit code 4: max. iterations reached
				infostruct.OutputMessage = outputMessage;
				infostruct = truncateInfoStruct(infostruct,i_iter);
				% Assign infostruct to second (optional) output
				varargout{2} = infostruct;
			end
		end
		% Game over: max. iterations reached
		return
	end
	% If an improvement has been obtained
	% with respect to the last point,
	if cost < cost_old
		% Assemble the gradient flag vector. This vector stores information
		% about the difference stencils to use during linearization:
		% -1 for backward difference,
		%  0 for central difference (default),
		% +1 for forward difference.
		% This way it is ensured that the evaluations during linearization
		% do not violate the bounds of x_u_min and x_u_max.
		gradient_flag = assembleGradientFlags(x_u,i_t_v,del_lin,x_u_min,x_u_max);
		% Linearize relevant subsystem at current operating point
		jaco = jj_lin (sys, x_u, n_x, i_t_v, i_t_r, del_lin, gradient_flag);
		% Assign current Jacobian to infostruct
		if nargout > 5
			infostruct.jaco{i_iter+1} = jaco;
		end
		% Singular Value Decomposition of the sensitivity matrix
		[u, s, v] = svd (jaco);
		% A singular value is assumed to be "zero", if it is 1e12 times smaller 
		% than the maximum singular value. Such a singular value indicates a rank deficiency.
		sv_min = s(1,1)*1e-12;
		% Find the indices of those "zero-singular-values"
		i_sv_zero = find (abs (diag (s)) <= sv_min);
		% If there are no zero-singular-values,
		if isempty(i_sv_zero)
			% accept and save this new point.
			% Important for a possible step size bisection later on
			x_u_old = x_u;
			% Save the cost value of this new point for a comparison later on
			cost_old = cost;
			% Reset step size bisection counter
			i_bisec = 0;
			% Assuming a linear system, the variation of the trim variables
			% necessary to compensate the trim requirements error can directly
			% be calculated by the inversion of the linear subsystem model
			% (differential equations and output equations)
			del_t_v = jaco\del_t_r;
			[del_t_v,stepLengthViolated] = forceConformanceWithMaximumStepLength(del_t_v,del_max(i_t_v));
			[del_t_v,hardBoundsViolated] = forceConformanceWithHardLimits(x_u,del_t_v,i_t_v,x_u_min,x_u_max);
			% Set the LimitedStep flag accordingly
			if stepLengthViolated && hardBoundsViolated
				LimitedStep = 3;
			elseif hardBoundsViolated
				LimitedStep = 2;
			elseif stepLengthViolated
				LimitedStep = 1;
			else
				LimitedStep = 0;
			end
			% If at least one element of the step is nonzero
			if ~all(del_t_v == 0)
				% Accept the step and calculate new trim point.
				% Always use old value *before* first bisection
				x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;
				% Disassemble the generalized input vector
				x_tr = x_u(1:n_x);
				u_tr = x_u(n_x+1:end);
				d_tr = d_y_tr(1:n_d);
				y_tr = d_y_tr(n_d+1:end);
				% Update step type: 1 = regular iteration
				if nargout > 5
					steptype = 1;
					infostruct.StepType(i_iter+1)       = steptype;
					infostruct.BisecCounter(1,i_iter+1) = i_bisec;
					infostruct.LimitedStep(1,i_iter+1)  = LimitedStep;
				end
			else
				% Game over: no step in the required direction possible (dead end)
				% Generate output message
				l1 = ['Trim failed in iteration ', num2str(i_iter),':'];
				l2 = 'No step in the required direction allowed because';
				l3 = 'the calculated step would violate the specified bounds.';
				l4 = ' ';
				l5 = 'Either the trim requirements cannot be met';
				l6 = 'or the specified min./max. bounds are too strict.';
				l7 = ' ';
				l8 = 'Try trimming without hard bounds';
				l9 = 'or try different initial values.';
				outputMessage = {l1; l2; l3; l4; l5; l6; l7; l8; l9};
				% Display error message if desired
				if EnableMessages
					errordlg (outputMessage,'Dead end');
				end
				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end
				if nargout > 4
					errflg = 1;
					varargout{1} = errflg;
					if nargout > 5
						infostruct.ExitCode = 5; % Exit code 5: dead-end
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);
						varargout{2} = infostruct;
						return
					end
				end
			end
		else % if the Jacobian is singular
			% Update step type and bisection counter
			if nargout > 5
				steptype = 3;
				infostruct.StepType(i_iter+1)       = steptype;
				infostruct.BisecCounter(1,i_iter+1) = i_bisec+1;
				infostruct.LimitedStep(1,i_iter+1)  = 0;
			end
			% If the initial Jacobian is the problem
			if i_iter == 0
				% Generate/display output message
				outputMessage = generateErrorMessage(x_nam,u_nam,d_nam,y_nam,u,v,i_t_v,i_t_r,sv_min,i_sv_zero,0,EnableMessages);
				% Update optional outputs and return
				if nargout > 4
					% Update error flag and return it as 5th parameter
					errflg = 1;
					varargout{1} = errflg;
					% If a 6th output parameter is used,
					if nargout > 5
						% Update infostruct and remove unassigned parts
						infostruct.ExitCode	= 1; % Exit code 1: singular initial Jacobian
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);
						% Store Jacobian (necessary because i_iter = 0)
						infostruct.jaco{1} = jaco;
						% Assign infostruct to second (optional) output
						varargout{2} = infostruct;
					end
				end
				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end
				% Game over: singular Jacobian at initial point
				return
			end
			% If step size has not been bisected N_bisec_max times yet,
			if i_bisec < N_bisec_max
				% bisect step size and change sign
				del_t_v = del_t_v/2;
				% and increment bisection counter
				i_bisec = i_bisec + 1;
				% Calculate new trim point.
				% Always use last value where Jacobian was regular and cost
				% had decreased
				x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;
			% If step size has already been bisected N_bisec_max times before,
			else
				% Prepare error message and abort
				% Assemble output message
				l1 = ['Maximum number of consecutive bisections (',num2str(N_bisec_max),') reached'];
				l2 = 'while trying to restore a regular Jacobian.';
				l3 = ['Program was aborted after ', int2str(i_iter), ' iteration(s)'];
				l4 = ['with a cost value of ', num2str(cost), '.'];
				l5 = 'Try different initial values.';
				l6 = 'Or try to reduce step sizes.';
				outputMessage = {l1; l2; l3; l4; l5; l6};
				% Display error message if desired
				if EnableMessages
					errordlg (outputMessage,'Program aborted');
				end
				% Release system
				if (CompileFlag == 1) || (CompileFlag == 2)
					feval (sys, [], [], [], 'term');
				end
				% Game over: max. number of bisections reached
				if nargout > 4
					errflg = 1;
					varargout{1} = errflg;
					if nargout > 5
						% Exit code 2: max. number of bisections reached
						% when trying to restore a regular Jacobian
						% Update infostruct and remove unassigned parts
						infostruct.ExitCode	= 2;
						infostruct.OutputMessage = outputMessage;
						infostruct = truncateInfoStruct(infostruct,i_iter);
						varargout{2} = infostruct;
					end
				end
				return
			end
		end
	% If the cost function has degraded
	else
		% Store step type 2 and bisection counter
		if nargout > 5
			steptype = 2;
			infostruct.StepType(i_iter+1)       = steptype;
			infostruct.BisecCounter(1,i_iter+1) = i_bisec+1;
			infostruct.LimitedStep(1,i_iter+1)  = 0;
			% Store Jacobian as empty matrix to make clear that it has not
			% been updated since last successful iteration step
			infostruct.jaco{i_iter+1} = [];
		end
		if i_bisec < N_bisec_max
			% bisect step size and change sign
			del_t_v = del_t_v/2;
			% and increment bisection counter
			i_bisec = i_bisec + 1;
			% Calculate new trim point.
			% Always use last value where Jacobian was regular and cost
			% had decreased
			x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;
		else
			% Prepare error message and abort
			% Assemble output message
			l1 = ['Maximum number of consecutive bisections (',num2str(N_bisec_max),') reached'];
			l2 = 'while trying to reduce the cost function.';
			l3 = ['Program was aborted after ', int2str(i_iter), ' iteration(s)'];
			l4 = ['with a cost value of ', num2str(cost), '.'];
			l5 = 'Try different initial values.';
			l6 = 'Or try to reduce step sizes.';
			outputMessage = {l1; l2; l3; l4; l5; l6};
			% Display error message if desired
			if EnableMessages
				errordlg (outputMessage,'Program aborted');
			end
			% Release system
			if (CompileFlag == 1) || (CompileFlag == 2)
				feval (sys, [], [], [], 'term');
			end
			% Game over: max. number of bisections reached
			if nargout > 4
				errflg = 1;
				varargout{1} = errflg;
				if nargout > 5
					% Exit code 3: max. number of bisections reached while
					% trying to reduce the cost function
					infostruct.ExitCode	= 3;
					infostruct.OutputMessage = outputMessage;
					infostruct = truncateInfoStruct(infostruct,i_iter);
					varargout{2} = infostruct;
				end
			end
			return
		end
	end
end
end
function [del_t_v,stepLengthViolated] = forceConformanceWithMaximumStepLength(del_t_v,del_t_v_max)
% This function ensures that the calculated step length del_t_v does not
% violate the maximum step length del_t_v_max.
	% Calculate maximum ratio between necessary and allowed trim
	% step sizes
	ratio_t_v = del_t_v ./ del_t_v_max;
	max_rat = max(abs(ratio_t_v));
	% If allowed step size has been exceeded,
	if max_rat > 1
		% scale all state and input step sizes, 
		% in order to exploit most of the allowed step size
		del_t_v = del_t_v/max_rat;
		stepLengthViolated = 1;
	else
		stepLengthViolated = 0;
	end
end
function [del_t_v,hardBoundsViolated] = forceConformanceWithHardLimits(x_u,del_t_v,i_t_v,x_u_min,x_u_max)
% This function ensures that the calculated step length del_t_v for the
% current x_u does not violate the bounds of x_u_min and x_u_max.
	hardBoundsViolated = 0;
	ratio_t_v = del_t_v./(x_u_max(i_t_v) - x_u(i_t_v));
	max_rat = max(ratio_t_v);
	if max_rat > 1
		del_t_v = del_t_v/max_rat;
		hardBoundsViolated = 1;
	end
    
	ratio_t_v = del_t_v./(x_u(i_t_v) - x_u_min(i_t_v));
	min_rat = min(ratio_t_v);
 
	if min_rat < -1
		del_t_v = del_t_v/-min_rat;
		hardBoundsViolated = 1;
	end
end
function gradient_flag = assembleGradientFlags(x_u,i_t_v,del_lin,x_u_min,x_u_max)
% This function assembles the vector which declares if the gradients of a
% trim variable should be calculated using central difference (flag = 0),
% forward difference (flag = 1) or backward difference (flag = -1). Central
% differences are used by default, and forward/backward differences are
% only used if the linearization step from the current working point would
% violate the upper or lower bounds x_u_max/x_u_min.
% Note: The implemented logic only works if all elements of the linearization
%       step vector del_lin are positive.
	n_tv = length(i_t_v);
	gradient_flag = zeros(n_tv,1);
	for ii = 1:n_tv
		if (x_u(i_t_v(ii)) + del_lin(i_t_v(ii))) > x_u_max(i_t_v(ii))
			gradient_flag(ii) = -1;
		elseif (x_u(i_t_v(ii)) - del_lin(i_t_v(ii))) < x_u_min(i_t_v(ii))
			gradient_flag(ii) = +1;
		else
			gradient_flag(ii) = 0;
		end
	end
end
function infostruct = truncateInfoStruct(infostruct_in,n_iter)
% This function truncates the fields of infostruct to the number of
% iterations actually performed. This is necessary since it was
% pre-allocated with a maximum size, but the final structure should
% contain only the actual iterations (n_iter plus the initial point).
	% Copy/paste structure - some entries will remain the same
	infostruct = infostruct_in;
	infostruct.n_iter = n_iter;
	infostruct.x = infostruct_in.x(:,1:n_iter+1);
	infostruct.u = infostruct_in.u(:,1:n_iter+1);
	infostruct.d = infostruct_in.d(:,1:n_iter+1);
	infostruct.y = infostruct_in.y(:,1:n_iter+1);
	infostruct.cost = infostruct_in.cost(1:n_iter+1);
	% Truncate step information
	infostruct.StepType		= infostruct_in.StepType(1:n_iter);
	infostruct.BisecCounter = infostruct_in.BisecCounter(1:n_iter);
	infostruct.LimitedStep  = infostruct_in.LimitedStep(1:n_iter);
	infostruct.jaco = infostruct_in.jaco(1:n_iter);
end
function errorMessage = generateErrorMessage(x_nam,u_nam,d_nam,y_nam,u,v,i_t_v,i_t_r,sv_min,i_sv_zero,bisectionWasUsed,enableMessages)
% This function generates the message that is displayed in case of a singular
% Jacobian matrix. It informs the user which trim variables and/or trim
% requirements might be the problem.
	% Generate top and bottom part of the error message to distinguish
	% whether or not bisection was used to restore a regular Jacobian
	% Note: the middle part of the message is common in both cases.
	if bisectionWasUsed == 0
		messageTitle = 'Singular initial Jacobian matrix!';
		upperPart = ' ';
		lowerPart = {...
			' ';
			'Choose different trim variables and/or trim requirements.';
			'Or try different initial values.';
			' ';
			'The returned trim point is not valid.';
			'If you are using trimmod, you can use the';
			'Untrim menu entry to return to the pre-trim state.'};
	else
		messageTitle = 'Singular Jacobian matrix!';
		upperPart = {...
			'Originally the trim problem was correctly';
			'set up (regular Jacobian matrix).';
			'This property was lost along the way and';
			'could not be restored by bisection.'};
		lowerPart = {...
			' ';
			'The returned trim point is not valid.';
			'If you are using trimmod, you can use the';
			'Untrim menu entry to return to the pre-trim state.'};
	end
	% Initialize the core part of the error message
	middlePart = [];
	% Assemble cell arrays containing the names of all trim variables and trim requirements
	trim_variables = [x_nam; u_nam];
	trim_requirements = [d_nam; y_nam];
	% Loop over all zero-singular-values
	for i_sv = i_sv_zero'
		% Find those elements of the corresponding singular vectors that are not "zero"
		u_sing = find (abs (u(:,i_sv)) > sv_min);
		v_sing = find (abs (v(:,i_sv)) > sv_min);
		% Separating empty line
		l10 = {' '};
		% If there is only one zero element in the left singular vector,
		if length (u_sing) == 1
			% prepare the corresponding error message
			l11 = {'The trim requirement'};
			l12 = trim_requirements(i_t_r(u_sing));
			l13 = {'could not be affected by any trim variable.'};
		% If there are more than one zero element in the left singular vector
		else
			% prepare the corresponding error message
			l11 = {'The trim requirements'};
			l12 = trim_requirements(i_t_r(u_sing));
			l13 = {'linearly depend on each other.'};
		end
		% Separating empty line
		l14 = {' '};
		% If there is only one zero element in the right singular vector,
		if length (v_sing) == 1
			% prepare the corresponding error message
			l15 = {'The trim variable'};
			l16 = trim_variables(i_t_v(v_sing));
			l17 = {'does not affect any trim requirement.'};
		% If there are more than one zero element in the right singular vector
		else
			% prepare the corresponding error message
			l15 = {'The trim variables'};
			l16 = trim_variables(i_t_v(v_sing));
			l17 = {'linearly depend on each other.'};
		end
		% Display error message
		if enableMessages
			errordlg([upperPart; l10; l11; l12; l13; l14; l15; l16; l17; lowerPart],messageTitle);
		end
		% Enhance the middle part of the complete error message
		middlePart = [middlePart; l10; l11; l12; l13; l14; l15; l16; l17;]; %#ok<AGROW>
	end
	% Assemble final output message
	errorMessage = [messageTitle; upperPart; middlePart; lowerPart];
end

function jaco = jj_lin (sys, x_u, n_x, i_x_u, i_d_y, del_x_u, gradient_flag)
%JJ_LIN   Subsystem linearisation of a nonlinear ordinary differential equation system 
% 
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y)  
%   linearizes the system with the name SYS at the operating point
%   which is defined by the generalized input vector X_U.
%   N_X is the length of the original state vector X. 
%   It is needed for the disassembling of the X_U vector 
%   in the parameter list of the system calls.
%   The Jacobian matrix JACO only contains the subsystem defined by the
%   index vectors I_X_U and I_D_Y.
%
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y, DEL_X_U)
%   additionally allows the specification of the perturbation levels
%   DEL_X_U to be used for the gradient calculations.
%
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y, DEL_X_U, GRADIENT_FLAG)
%   additionally allows the specification of the difference stencil to use.
%   GRADIENT_FLAG can be defined as a vector of N = length(i_x_u) elements.
%   The values must be {-1,0,1}:
%   -1: backward differences,
%    0: central differences (default),
%   +1: forward differences.
%   GRADIENT_FLAG can also be passed as a scalar, in which case the same
%   difference stencil is applied for all trim variables.
%
%   Function dependencies: getGeneralizedOutputVector (for model evaluation)
%
%   Copyright 2000-2004, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de
%   Version 1.3, 22.07.2020
%   Daniel Kiehn, German Aerospace Center (DLR), daniel.kiehn@dlr.de
%   - added pre-allocation of the Jacobian matrix
%   - index vectors can now be column or row vectors
%   - added optional input gradient_flag
%   - Moved the model evaluation to new function getGeneralizedOutputVector
% If the user did not define the perturbation levels
if nargin < 6
	% Use default perturbation levels
	del_x_u = 1e-6*(1 + abs(x_u)); 
end
% If the user did not provide the gradient flag
if nargin < 7
	gradient_flag = zeros(size(i_x_u)); % Use central differences by default
else
	if isscalar(gradient_flag)
		gradient_flag = gradient_flag*ones(size(i_x_u));
	elseif length(gradient_flag) ~= length(i_x_u)
		errmsg = ['The length of gradient_flag (',num2str(length(gradient_flag)),')',...
				  ' does not match the length of i_x_u (',num2str(length(i_x_u)),').'];
		error(errmsg);
	end
end
% Determine vector lengths
n_i_x_u = length (i_x_u);
n_i_d_y = length (i_d_y);
% Set time to zero explicitly (assume a time-invariant system)
t = 0;
% Pre-allocate Jacobian matrix
jaco = zeros(n_i_d_y,n_i_x_u);
% Loop over all generalized inputs: generate one column of the Jacobian for
% every generalized input to be linearized and build up the matrix columnwise
for ii = 1:n_i_x_u
	% Save whole generalized input vector in x_u_left and x_u_right,
	% because we only want to waggle some specific generalized inputs
	x_u_left  = x_u;
	x_u_right = x_u;
	% Waggle the ii-th generalized input
	x_u_left(i_x_u(ii))  = x_u(i_x_u(ii)) - del_x_u(i_x_u(ii));
	x_u_right(i_x_u(ii)) = x_u(i_x_u(ii)) + del_x_u(i_x_u(ii));
	% Central difference
	if gradient_flag(ii) == 0
		d_y_left  = getGeneralizedOutputVector(sys, x_u_left,  n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u_right, n_x, t);
		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(2*del_x_u(i_x_u(ii)));
	% Right-sided difference
	elseif gradient_flag(ii) == 1
		d_y_left  = getGeneralizedOutputVector(sys, x_u,       n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u_right, n_x, t);
		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(del_x_u(i_x_u(ii)));
	% Left-sided difference
	elseif gradient_flag(ii) == -1
		d_y_left  = getGeneralizedOutputVector(sys, x_u_left, n_x, t);
		d_y_right = getGeneralizedOutputVector(sys, x_u,      n_x, t);
		jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(del_x_u(i_x_u(ii)));
	else
		errmsg = ['Wrong value of gradient_flag for generalized input #',...
				  num2str(ii),': ',num2str(gradient_flag(ii)),'. ',...
				  'Only {-1,0,1} are accepted!'];
		error(errmsg);
	end
	jaco(:,ii) = jaco_column;
end
end

function d_y = getGeneralizedOutputVector(sys,x_u,n_x,t)
% function d_y = getGeneralizedOutputVector(sys,x_u,n_x,t)
%
% Evaluates the generalized output vector (state derivatives and outputs)
% of the Simulink model sys at the working point x_u at time t.
%
% Inputs:
% sys       Name of the Simulink model (without file extension) as string
% x_u		Generalized input vector: column vector of states x followed by
%           inputs u
% n_x       Number of states. This is needed to disassemble x_u into x and u
% t         Time in seconds
%
% Outputs:
% d_y		Generalized output vector: column vector of state derivatives d
%           followed by outputs y
%
% 22.07.2020, Daniel Kiehn, German Aerospace Center (DLR), daniel.kiehn@dlr.de
% Calculate outputs and derivatives at the working point (x,u).
% Important: We have to calculate the outputs first!
%            The derivatives would be wrong otherwise.
% Unbelievable but true: Mathworks argues that this is not a bug but a feature!
% And furthermore: Sometimes it is even necessary to do the output calculation twice
% before you get the correct derivatives!
% Mathworks says, they will take care of that problem in one of the
% next releases...
y = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 'outputs');
y = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 'outputs');
d = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 'derivs');
d_y = [d; y];
end
