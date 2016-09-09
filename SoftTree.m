classdef SoftTree < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        X;
        Y;
        V;
        R;
        type;
        treeRoot;
    end
    
    methods
        function this = SoftTree(X,Y,V,R)
            this.treeRoot = Node();            
            this.treeRoot.w = zeros(1:length(X(1,:)));
            this.treeRoot.w0 = 0;
            
            this.X = X;
            this.Y = Y;
            
            this.V = V;
            this.R = R;            
            for i = 1:length(Y)
                this.treeRoot.w0 = this.treeRoot.w0 +Y(i);
            end
            this.treeRoot.w0 = this.treeRoot.w0/length(Y);
            
            this.type = 'c';
        end
        
        function this = train(this)            
            this.treeRoot.splitNode(this, this.X, this.Y, this.V, this.R); % THISSS
        end
        
        function ret = evaluate(this, x)
            if (this.type == 'r')
                ret = treeRoot.evaluate(x);
                return;
            else (this.type == 'c');
                lamp = this.treeRoot.evaluate(x);
                ret = cppsigmoid(lamp);
                return;
            end            
        end
        
        function ret = meanSqErr(X,Y)
            err = 0;
            for i = 1:length(Y)
                y = this.evaluate(X(i,:));
                err = err + (Y(i)-y)*(Y(i)-y);
            end
            err = err/length(Y);
            ret = err;
            return;
        end
        function ret = errRate(this, X,Y)
            err = 0;
            for i = 1:length(Y)
                y = this.evaluate(X(i,:)); 
                err = err + (Y(i) ~= (y > 0.5));
            end
            ret = err/length(Y);
            return;
        end
        
        function ret = size(this)
            ret = this.treeRoot.size();
            return;
        end
        
        function ret = print(this)
            ret = this.treeRoot.print(1);
            return;
        end
        
        
    end
    
end

