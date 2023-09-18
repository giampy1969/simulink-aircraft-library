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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jj_3m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ...
  [x_tr, u_tr, d_tr, y_tr] = ...
  jj_3m (...
  sys, ...
  x, u, d, y, ...
  i_x, i_u, i_d, i_y, ...
  x_nam, u_nam, d_nam, y_nam, ...
  del_x_max, del_u_max, ...
  del_x_lin, del_u_lin, ...
  options)

%JJ_TRIM   Trim point determination of a nonlinear ordinary differential equation system
%
%   [X_TR, U_TR, D_TR, Y_TR] = JJ_TRIM (SYS, X, U, D, Y, I_X, I_U, I_D, I_Y, X_NAM, U_NAM, D_NAM, Y_NAM)
%   trims the system SYS towards an operating point defined by the elements of
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
%   To see more help, enter TYPE JJ_TRIM.

%   [X_TR, ...] = JJ_TRIM (..., Y_NAM, DEL_X_MAX, DEL_U_MAX)
%   allows the additional specification of maximum alterations
%   of state and input during one trim step.
%   The lengths of DEL_X_MAX and DEL_U_MAX are equal to those of
%   X and U respectively.
%   The default values of DEL_X_MAX and DEL_U_MAX are 1e42.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_MAX, DEL_X_LIN, DEL_U_LIN)
%   allows the additional specification of the state and input step size
%   to be used for the calculation of the Jacobian-matrix (sensitivity matrix)
%   in the linearization procedure.
%   The default values of DEL_X_LIN and DEL_U_LIN are
%   1e-6*(1 + abs(x)) and 1e-6*(1 + abs(u)) respectively.
%
%   [X_TR, ...] = JJ_TRIM (..., DEL_U_LIN, OPTIONS)
%   allows the additional specification of
%   the maximum number of iterations "OPTIONS(1)" (default: 42) and
%   the cost value "OPTIONS(2)" (default: 1e-9) to be gained.

%   Copyright 2000-2004, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.2     26.05.2000
%
%   The names of inputs, ... outputs for the error messages
%   are now transferred via the parameter list.
%
%   The precompilation and the release of the system is now done in JJ_TRIM.


% Feed through all initial vectors, 
% usefull in case of an emergency exit
x_tr = x;
u_tr = u;
d_tr = d;
y_tr = y;

% Determine lengths of basic vectors
n_x  = length (x);
n_u  = length (u); 
n_d  = length (d);
n_y  = length (y);

n_i_x = length (i_x);
n_i_u = length (i_u);
n_i_y = length (i_y);
n_i_d = length (i_d);

% Assemble generalized input vector and generalized output vector
x_u = [x; u];
d_y = [d; y];

% Determine length of generalized input vector and generalized output vector
n_x_u = length (x_u);
n_d_y = length (d_y);

% Assemble trim variable index vector and trim requirement index vector
i_t_v = [i_x; i_u + n_x]; 
i_t_r = [i_d; i_y + n_d];

% Determine number of trim variables and trim requirements
n_t_v = n_i_x + n_i_u;
n_t_r = n_i_y + n_i_d;

% There have to be as many trim variables as there are trim requirements
if n_t_r ~= n_t_v
  
  l1 = ['The number of trim variables: ', int2str(n_t_v)];
  l2 = 'does not equal';
  l3 = ['the number of trim requirements: ', int2str(n_t_r)];
  l4 = ' ';
  l5 = 'The returned trim point is not valid.';
  h1 = 'Error';
  %  errordlg ({l1; l2; l3; l4; l5}, h1); 
  disp([h1,': ',l1,' ',l2,' ',l3,' ',l4,' ',l5]);
  
  % Game over
  return
  
end

% There should be at least one trim variable and one trim requirement
if ~n_t_r
  
  l1 = 'There should be at least';
  l2 = 'one trim variable and';
  l3 = 'one trim requirement.';
  h1 = 'Nothing to trim';
  
  disp([h1,': ',l1,' ',l2,' ',l3]);
  % warndlg ({l1; l2; l3}, h1);
  
end

% If no maximum step sizes have been defined,
if nargin < 14
  
  % set defaults
  del_max = 1.e42*ones (n_x_u, 1); 
  
  % otherwise assemble generalized maximum step vector
else
  
  del_max = [del_x_max; del_u_max];
  
end

% If no step sizes for the linearization have been defined,
if nargin < 16
  
  % set defaults
  del_lin = 1e-6*(1 + abs(x_u)); 
  
  % otherwise assemble generalized linearization step vector
else
  
  del_lin = [del_x_lin; del_u_lin];
  
end

% If no options have been defined, 
if nargin < 18
  
  % set defaults (maximum number of iterations and cost value to be gained)
  options(1) = 42;
  options(2) = 1e-9;
  
end

% Save and rename the options 
n_iter   = options(1);
cost_tbg = options(2);

% Set time to zero explicitely (assume a time invariant system)
t = 0;

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
feval (sys, [], [], [], 'compile');
    
% Loop over maximum n_iter iterations
for i_iter = 0 : n_iter
  
  % Calculate outputs and derivatives at the current trim point.
  % Important: We have to calculate the outputs first!
  %            The derivatives would be wrong otherwise.
  % Unbelievable but true: Mathworks argues that this is not a bug but a feature!
  % And furthermore: Sometimes it is even necessary to do the output calculation twice
  % before you get the correct derivatives!
  % Mathworks says, they will take care of that problem in one of the next releases...
  y_tr = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 3); 
  y_tr = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 3); 
  d_tr = feval (sys, t, x_u(1:n_x), x_u(n_x+1:end), 1);
  d_y_tr = [d_tr; y_tr];
  
  % Calculate differences between required and current generalized output vectors
  del_d_y = d_y - d_y_tr;
  
  % Pick trim requirements only
  del_t_r = del_d_y(i_t_r);
  
  % Cost value is the maximum element of the trim requirement error vector
  cost = max (abs (del_t_r));
  
  % Cost is an empty matrix if there are no trim variables and trim requirements
  if isempty (cost)
    cost = 0;
  end
  
  % If current cost value has become smaller 
  % than the cost value to be gained
  if cost < cost_tbg
    
    % Output cost value and number of iterations needed
    l1 = ['A cost value of ',num2str(cost)];
    l2 = ['has been gained after ', int2str(i_iter), ' iteration(s).'];   
    h1 = 'Success';
    disp([h1,': ',l1,' ',l2]);
    % msgbox ({l1, l2}, h1);
    
    % Release system
    feval (sys, [], [], [], 'term');
    
    % Game over
    return
    
  end
  
  % If an improvement has been obtained
  % with respect to the last point,
  if cost < cost_old
    
    % accept and save this new point.
    % Important for a possible step size bisection later on
    x_u_old = x_u;
    
    % Save the cost value of this new point for a comparison later on
    cost_old = cost;    
    
    % Reset step size bisection counter
    i_bisec = 0;
    
    % Linearize relevant subsystem at current operating point
    jaco = jj_lin (sys, x_u, n_x, i_t_v, i_t_r, del_lin);
    
    % Singular Value Decomposition of the sensitivity matrix
    [u, s, v] = svd (jaco);
    
    % A singular value is assumed to be "zero", if it is 1e12 times smaller 
    % than the maximum singular value. Such a singular value indicates a rank deficiency.
    sv_min = s(1,1)*1e-12;
    
    % Find the indices of those "zero-singular-values"
    i_sv_zero = find (abs (diag (s)) <= sv_min);
    
    % If there are any zero-singular-values,
    if ~isempty (i_sv_zero) 
      
      % the jacobian matrix is singular.
      h1 = 'Singular Jacobian-Matrix';
      
      % Assemble cell arrays containing the names of all trim variables and trim requirements
      trim_variables = [x_nam; u_nam];
      trim_requirements = [d_nam; y_nam];
      
      % Loop over all zero-singular-values
      for i_sv = i_sv_zero'
        
        % Find those elements of the corresponding singular vectors that are not "zero"
        u_sing = find (abs (u(:,i_sv)) > sv_min);
        v_sing = find (abs (v(:,i_sv)) > sv_min);   
        
	      % Separating empty line
    	  l0 = {' '};

        % If there is only one zero element in the left singular vector,
        if length (u_sing) == 1
          
          % prepare the corresponding error message
          l1 = {'The trim requirement'};
          l2 = trim_requirements(i_t_r(u_sing));
          l3 = {'could not be affected by any trim variable.'};       
          
        % If there are more than one zero element in the left singular vector
        else
          
          % prepare the corresponding error message
          l1 = {'The trim requirements'};
          l2 = trim_requirements(i_t_r(u_sing));
          l3 = {'linearly depend on each other.'};       
          
        end 
        
        % Separating empty line
        l4 = {' '};
        
        % If there is only one zero element in the right singular vector,
        if length (v_sing) == 1
          
          % prepare the corresponding error message
          l5 = {'The trim variable'};
          l6 = trim_variables(i_t_v(v_sing));
          l7 = {'does not affect any trim requirement.'};       
          
        % If there are more than one zero element in the right singular vector
        else
          
          % prepare the corresponding error message
          l5 = {'The trim variables'};
          l6 = trim_variables(i_t_v(v_sing));
          l7 = {'linearly depend on each other.'};       
          
        end 
        
        l8 = {'Chose different trim variables and/or trim requirements.'};
        l9 = {'Or try different initial values.'};
        
        % Separating empty line
        l10 = {' '};
        
        l11 = {'The returned trim point is not valid.'};
        l12 = {'You can use the Untrim menu entry to return to the pre-trim state.'};
        
        % Output error message
        % errordlg ([l0; l1; l2; l3; l4; l5; l6; l7; l4; l8; l9; l10; l11; l12], h1);
        disp([h1,': ',l0{:},' ',l1{:},' ',l2{:},' ',l3{:},' ',l4{:},' ',l5{:},' ',l6{:},' ',l7{:},' ',l8{:},' ',l9{:}]);
        
      end
      
      % Release system
      feval (sys, [], [], [], 'term');
    
      % Game over
      return
      
    end
    
    % Assuming a linear system, the alteration of the trim variables 
    % necessary to compensate the trim requirements error can directly 
    % be calculated by the inversion of the linear subsystem model
    % (differential equations and output equations)
    del_t_v = jaco\del_t_r;
    
    % Calculate maximum ratio between allowed and necessary trim step size
    ratio_t_v = del_t_v ./ del_max(i_t_v);
    max_rat = max (abs (ratio_t_v));
    
    % If allowed step size has been exceeded,
    if max_rat > 1
      
      % scale all state and input step sizes, 
      % in order to exploit most of the allowed step size
      del_t_v = del_t_v/max_rat;
      
    end
    
    % If no improvement has been obtained
    % with respect to the last point,
  else
    
    % and if step size has not been bisected 200 times before,
    if i_bisec < 200
      
      % bisect step size and change sign
      del_t_v = -del_t_v/2;
      
      % and increment bisection counter
      i_bisec = i_bisec + 1;
      
      % If step size has already been bisected 200 times before,
    else
      
      % output error message and stop program
      l1 = 'Step size has been bisected 200 times.';
      l2 = ['Program was aborted after ', int2str(i_iter), ' iteration(s)'];   
      l3 = ['with a cost value of ', num2str(cost)];
      l4 = 'Try different initial values.';
      h1 = 'Program aborted';
      % errordlg ({l1; l2; l3; l4}, h1);
      disp([h1,': ',l1,' ',l2,' ',l3,' ',l4]);
      
      % Release system
      feval (sys, [], [], [], 'term');
    
      % Game over
      return
      
    end
    
  end
  
  % Calculate new trim point.
  % Always use old value *before* first bisection
  x_u(i_t_v) = x_u_old(i_t_v) + del_t_v;
  
  % Disassemble the generalized input vector
  x_tr = x_u(1:n_x);
  u_tr = x_u(n_x+1:end);
  
end

% If maximum number of iterations has been exceeded,
% output error message and abort program
l1 = ['Maximum number of iterations exceeded: ', int2str(i_iter)];   
l2 = 'Program was aborted';   
l3 = ['with a cost value of: ', num2str(cost)];
h1 = 'Program aborted';
% errordlg ({l1; l2; l3}, h1);
disp([h1,': ',l1,', ',l2,' ',l3]);

% Release system
feval (sys, [], [], [], 'term');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% jj_3m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ...
  jaco = ...
  jj_lin (...
  sys, ...
  x_u, n_x, ...
  i_x_u, i_d_y, ...
  del_x_u)

%JJ_LIN   Subsystem linearisation of a nonlinear ordinary differential equation system 
% 
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y)  
%   linearizes the system with the name SYS at the operating point,
%   that is defined by the generalized input vector X_U.
%   N_X is the length of the original state vector X. 
%   It is needed for the dissassembling of the X_U vector 
%   in the parameter list of the system calls.
%   The matrix JACO only contains the subsystem defined by the 
%   index vectors I_X_U and I_D_Y.
%
%   JACO = JJ_LIN (SYS, X_U, N_X, I_X_U, I_D_Y, DEL_X_U)
%   additionally allows the specification of the perturbation levels 
%   DEL_X_U to be used for the gradient calculations.

%   Copyright 2000-2004, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.2     26.05.2000


% If the user did not define the perturbation levels
if nargin < 6
  
  % Use default perturbation levels
  del_x_u = 1e-6*(1 + abs(x_u)); 
  
end

% Determine vector lengths
n_i_x_u = length (i_x_u);
n_i_d_y = length (i_d_y);

% Set time to zero explicitely (assume a time invariant system)
t = 0;

% Initialize matrices (will be constructed columnwise later on) 
jaco = [];

% Loop over all generalized inputs to be linearized.
% IMPORTANT: Vector has to be a row vector!
for i = i_x_u'
  
  % Save whole generalized input vector in x_u_left and x_u_right,
  % because we only want to waggle some specific generalized inputs 
  x_u_left  = x_u;
  x_u_right = x_u;
  
  % Waggle one generalized input
  x_u_left(i)  = x_u(i) - del_x_u(i);
  x_u_right(i) = x_u(i) + del_x_u(i);
  
  % Calculate outputs and derivatives at the current trim point.
  % Important: We have to calculate the outputs first!
  %            The derivatives would be wrong otherwise.
  % Unbelievable but true: Mathworks argues that this is not a bug but a feature!
  % And furthermore: Sometimes it is even necessary to do the output calculation twice
  % before you get the correct derivatives!
  % Mathworks says, they will take care of that problem in one of the next releases...
  y_left = feval (sys, t, x_u_left(1:n_x), x_u_left(n_x+1:end), 3);
  y_left = feval (sys, t, x_u_left(1:n_x), x_u_left(n_x+1:end), 3);
  d_left = feval (sys, t, x_u_left(1:n_x), x_u_left(n_x+1:end), 1);
  
  y_right = feval (sys, t, x_u_right(1:n_x), x_u_right(n_x+1:end), 3);
  y_right = feval (sys, t, x_u_right(1:n_x), x_u_right(n_x+1:end), 3);
  d_right = feval (sys, t, x_u_right(1:n_x), x_u_right(n_x+1:end), 1);
  
  % Assemble generalized output vectors
  d_y_left = [d_left; y_left];
  d_y_right = [d_right; y_right];
  
  % Generate one column of the jacobi-matrix for every generalized input to be linearized
  % and build up the matrix columnwise
  jaco_column = (d_y_right(i_d_y) - d_y_left(i_d_y))/(2*del_x_u(i));   
  jaco = [jaco, jaco_column]; 
  
end

% Copyright 2018 The MathWorks, Inc.