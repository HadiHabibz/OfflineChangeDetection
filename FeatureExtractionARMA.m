classdef FeatureExtractionARMA < FeatureExtractionMethod
    
    properties(Access = protected)
        orderAR;
        orderMA;
    end % class datamembers 
    
    methods(Access = public)
        
        function obj = FeatureExtractionARMA()
            % Class constructor
            obj@FeatureExtractionMethod();            
        end % function FeatureExtractionMisc
        
        function order = getOrderAR(obj)
            order = obj.orderAR;
        end
        
        function order = getOrderMA(obj)
            order = obj.orderMA;
        end 
                           
    end % class public services
    
  methods(Access = protected)          
            
      function testvalues = computeFeatures(obj, Y)  
          
          % Do nothing!
          testvalues = Y;
            
      end % class computeFeatures 
                               
        function obj = performInitializations(obj)
            % Empty
            
        end % function performOtherPrecomputations     
        
        function obj = parseInputVariables(obj, options)

            obj = parseInputVariables@FeatureExtractionMethod(...
                obj, options); 
            
            parser = inputParser;
            parser.KeepUnmatched = true;
            addParameter(parser, 'orderAR', 1);
            addParameter(parser, 'orderMA', 0);
            parse(parser, options{:});

            obj = setOrderAR(obj, parser.Results.orderAR);
            obj = setOrderMA(obj, parser.Results.orderMA);

        end % function parseInputVariables
        
        function obj = setOrderAR(obj, order)
            % Setter function for the order of AR model
            
            if ~isnumeric(order) || length(order) ~= 1
                error('Model orders must be a natural number. ');
            end
            
            if order <= 0
                error('Model orders must be a natural number. ');
            end
            
            obj.orderAR = order;
            
        end % function setOrder
        
        function obj = setOrderMA(obj, order)
            
            if ~isnumeric(order) || length(order) ~= 1
                error('Model orders must be a non-negative integer. ');
            end
            
            if order < 0
                error('Model orders must be a non-negative integer. ');
            end
            
            obj.orderMA = order;
            
        end % function setOrderMA
        
    end % class utilities
    
    
end % class featureExtractionARMA