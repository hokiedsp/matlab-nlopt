classdef options < matlab.mixin.Copyable & matlab.mixin.SetGet & matlab.mixin.CustomDisplay
   %nlopt.options   NLopt option class
   %
   %
   %   See also
   
   properties (Dependent,SetAccess=private)
      Algorithm
      Dimension
   end
   
   properties (Dependent)
      ObjectiveLimit   % StopVal
      FunctionRelativeTolerance % FTolRel
      StepRelativeTolerance % XTolRel
      FunctionTolerance % FTolAbs
      StepTolerance % XTolAbs
      MaxFunctionEvaluations % MaxFunEvals
      MaxEvaluationDuration  %MaxEvalTime
      Population
      VectorStorage
      InitialStepSize
   end
   
   properties
      HessMultFcn % W = HessMultFcn(H,v);
%       Display %
      OutputFun % function_handle withs signature: stop = outfun(x,optimValues,state)
      
      SubproblemAlgorithm  % scalar nlopt.options for local optimization configuration (used only by some global optimization routines)
      
      ConstraintTolerance = 1e-12;
   end
   
   properties (Access = protected, Hidden, NonCopyable, Transient)
      backend % Handle to the backend C++ class instance
   end
   
   methods (Access = protected, Static, Hidden)
      varargout = mexfcn(varargin) % mexCopyableObjectHandler mex function, defined by derived class
   end
   
   methods
      function obj = options(algorithm, dim, varargin)
         validateattributes(algorithm,{'char'},{'row'},mfilename,'algorithm');
         validateattributes(dim,{'double'},{'scalar','integer','positive'},mfilename,'dim');
         
         % instantiate mex backend
         obj.mexfcn(obj, algorithm, dim);

         % append default stopping condition
         if nargin>2
            if strcmpi(varargin{1},'subproblem')
               varargin(1) = [];
            else
               varargin = [{'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
                  'MaxFunctionEvaluations',100*dim},varargin];
            end
         end
         
         % set properties
         if ~isempty(varargin)
            obj.set(varargin{:});
         end
      end
      
      function delete(obj)
         % destroy the mex backend
         if ~isempty(obj.backend)
            obj.mexfcn(obj, 'delete');
         end
      end
   end
   
   methods (Static)
      function ver = getNLoptVersion()
         ver = nlopt.options.mexfcn('getNLoptVersion');
      end
      function varargout = getAlgorithms()
         names = nlopt.options.mexfcn('getAlgorithms',nargout==0);
         names(:) = sortrows(names);
         if nargout > 0
            varargout{1} = names;
         else
            tbl = cell2table(names,'VariableNames',{'Name','Description'});
            disp(tbl);
            disp('For the details of these algorithm, see <a href="https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/">the official NLopt Documentation</a>');
         end
      end
   end
   
   methods(Access = protected)
      % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all four properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         
         % Make a deep copy of the C++ object
         obj.mexfcn(obj, 'copy', cpObj);
      end
   end
   
   methods
      function val = get.Algorithm(obj)
         val = obj.mexfcn(obj,'getAlgorithm');
      end
      function val = get.Dimension(obj)
         val = obj.mexfcn(obj,'getDimension');
      end
      function val = get.ObjectiveLimit(obj)   % StopVal
         val = obj.mexfcn(obj,'getFunctionStopValue');
      end
      function set.ObjectiveLimit(obj,val)
         validateattributes(val,{'double'},{'scalar','nonnan'});
         obj.mexfcn(obj,'setFunctionStopValue',val);
      end
      function val = get.FunctionRelativeTolerance(obj) % FTolRel
         val = obj.mexfcn(obj,'getFunctionRelativeTolerance');
         if val<=0
            val = 'off';
         end
      end
      function set.FunctionRelativeTolerance(obj,val)
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','finite'});
         end
         obj.mexfcn(obj,'setFunctionRelativeTolerance',val);
      end
      function val = get.StepRelativeTolerance(obj) % FTolRel
         val = obj.mexfcn(obj,'getStepRelativeTolerance');
         if val<=0
            val = 'off';
         end
      end
      function set.StepRelativeTolerance(obj,val) % XTolRel
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','finite'});
         end
         obj.mexfcn(obj,'setStepRelativeTolerance',val);
      end
      function val = get.FunctionTolerance(obj) % FTolAbs
         val = obj.mexfcn(obj,'getFunctionAbsoluteTolerance');
         if val<=0
            val = 'off';
         end
      end
      function set.FunctionTolerance(obj,val)
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','finite'});
         end
         obj.mexfcn(obj,'setFunctionAbsoluteTolerance',val);
      end
      function val = get.StepTolerance(obj) % XTolAbs
         val = obj.mexfcn(obj,'getStepAbsoluteTolerance');
         if isscalar(val) && val<=0
            val = 'off';
         end
      end
      function set.StepTolerance(obj,val)
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            if isscalar(val)
               validateattributes(val,{'double'},{'scalar','positive','finite'});
            else
               validateattributes(val,{'double'},{'vector','numel',obj.Dimension,'nonnegative','finite'});
            end
         end
         obj.mexfcn(obj,'setStepAbsoluteTolerance',val);
      end
      function val = get.MaxFunctionEvaluations(obj) % MaxEval
         val = obj.mexfcn(obj,'getMaxFunctionEvaluations');
         if val<=0
            val = 'off';
         end
      end
      function set.MaxFunctionEvaluations(obj,val)
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','finite'});
         end
         obj.mexfcn(obj,'setMaxFunctionEvaluations',val);
      end
      function val = get.MaxEvaluationDuration(obj) % MaxTime
         val = obj.mexfcn(obj,'getMaxEvaluationDuration');
         if val<=0
            val = 'off';
         end
      end
      function set.MaxEvaluationDuration(obj,val)
         try
            validatestring(val,{'off'});
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','finite'});
         end
         obj.mexfcn(obj,'setMaxEvaluationDuration',val);
      end
      function val = get.Population(obj)
         val = obj.mexfcn(obj,'getMaxEvaluationDuration');
         if val==0
            val = 'auto';
         end
      end
      function set.Population(obj,val)
         try
            validatestring(val,{'auto'})
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','integer'});
         end
         obj.mexfcn(obj,'setVectorStorage',val);
      end
      function val = get.VectorStorage(obj)
         val = obj.mexfcn(obj,'getVectorStorage');
         if val==0
            val = 'auto';
         end
      end
      function set.VectorStorage(obj,val)
         try
            validatestring(val,{'auto'})
            val = 0;
         catch
            validateattributes(val,{'double'},{'scalar','positive','integer'});
         end
         obj.mexfcn(obj,'setVectorStorage',val);
      end
      function val = get.InitialStepSize(obj)
         val = obj.mexfcn(obj,'getInitialStepSize');
      end
      function set.InitialStepSize(obj,val)
         validateattributes(val,{'double'},{'scalar','positive','finite'});
         obj.mexfcn(obj,'setInitialStepSize',val);
      end
      
      function set.HessMultFcn(obj,val)
         validateattributes(val,{'function_handle'},{'scalar'});
         obj.HessMultFcn = val;
      end
      
      function set.OutputFun(obj,val)
         % stop = outfun(x,optimValues,state);
         validateattributes(val,{'function_handle'},{'scalar'});
         obj.OutputFun = val;
      end
      
      function set.SubproblemAlgorithm(obj,val)
         validateattributes(val,{'nlopt.options'},{'scalar'});
         obj.SubproblemAlgorithm = val;
      end

      function set.ConstraintTolerance(obj,val)
         validateattributes(val,{'double'},{'scalar','nonnegative','finite'});
         obj.ConstraintTolerance = val;
      end
   end
   
   methods (Hidden)
      function varargout = fminunc(obj,fun,x0)
         % to be called only by nlopt.fminunc
         [varargout{1:nargout}] = obj.mexfcn(obj,'fminunc',fun,x0);
         if nargout>3
            varargout{4} = struct('funccount',varargout{4},'Algorithm',obj.Algorithm);
         end
      end
      
      function varargout = fmincon(obj,varargin)
         % to be called only by nlopt.fmincon
         [varargout{1:nargout}] = obj.mexfcn(obj,'fmincon',varargin{:});
         if nargout>3
            varargout{4} = struct('funccount',varargout{4},'Algorithm',obj.Algorithm);
         end
      end
   end
end
