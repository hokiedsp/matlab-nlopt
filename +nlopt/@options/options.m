classdef (Abstract) options < matlab.mixin.Copyable
%nlopt.options   nlopt options class
%
%
%   See also 
   
   properties (Access = protected, Hidden, NonCopyable, Transient)
      backend % Handle to the backend C++ class instance
   end
   
   methods (Access = protected, Static, Hidden, Abstract)
      varargout = mexfcn(varargin) % mexCopyableObjectHandler mex function, defined by derived class
   end
   
   methods
      function obj = options(varargin)
         % instantiate mex backend
         obj.mexfcn(obj, varargin{:});
      end
      
      function delete(obj)
         % destroy the mex backend
         if ~isempty(obj.backend)
            obj.mexfcn(obj, 'delete');
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
end
