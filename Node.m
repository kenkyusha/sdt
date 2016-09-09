    classdef Node < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        parent;
        left;
        right;
        isLeaf;
        isLeft;
        w;
        w0;
        y;
        v;
    end
    methods (Access = public)
        %
        % Constructor
        %
        function this = Node()            
            this.isLeaf = true;
            this.parent = 0;
            this.left = 0;
            this.right = 0;
            this.isLeft = true;
        end
        
    
        function ret = evaluate(this, x)
            if (this.isLeaf == true)
                this.y = this.w0; 
                ret = this.y;
                return;
            else
                this.v = cppsigmoid(cppdot(this.w,x)+this.w0);
                dummyLeft = this.left.evaluate(x);
                dummyLeft = this.v * dummyLeft;
                dummyRight = this.right.evaluate(x);
                dummyRight = dummyRight * (1-this.v);
                this.y = dummyLeft + dummyRight;
%                 this.y = this.v *(this.left.evaluate(x)) + (1-this.v)*(this.right.evaluate(x));
                ret = this.y;
                return;
            end 
        end
        
        function size = size(this)
            if (this.isLeaf == true)
                size = 1;
            else
               size =  (1 + this.left.size() + this.right.size());
            end
        end
        function print(this, depth)
            for i = 1:depth-1
                display('__');
            end
            if (this.isLeaf == false)
                for i = 1:length(this.w)
                    o1 = [this.w(i), ','];
                    disp(o1);
                end
            end
            disp(this.w0);
            if (this.isLeaf == false)
                this.left.print(depth+1);
                this.right.print(depth+1);
            end
        end
        
        function hardInit(this, X, Y)
            sv = [];
%             while (isempty(this.w) == 0) %kind of assertion
                total = 0;
                for i = 1:length(X);
                    t = 1;
                    m = this; 
                    p = Node();
                    while m.parent ~= 0
                        p = m.parent;
                        if (m.isLeft)
                            t = t* (cppsigmoid(cppdot(p.w, X(i,:)) + p.w0));
                        else
                            t = t* (1-(cppsigmoid(cppdot(p.w, X(i,:)) + p.w0)));
                        end
                        m = m.parent;
                    end
                    sv(i) = t;
                    total = total+t;
                end
                if (total <= 1)
                    this.w = zeros(1,size(X,2));
                    for i = 1:length(this.w)
                        this.w(i) = 0.005+(-0.005-0.005).*rand(1,1); %not sure whether this is correct
                    end
                    this.w0 = 0.005+(-0.005-0.005).*rand(1,1);
                    this.left.w0 = 0.005+(-0.005-0.005).*rand(1,1);
                    this.right.w0 = 0.005+(-0.005-0.005).*rand(1,1);
                end
                %dim = pos numbers only
                bestDim = -1;
                errBest = -1;

                numFeat = size(X,2);

                %look for best hardsplit
                for dim = 1:numFeat
                    f = [];
                    for i = 1:length(X);
                        f(i,1) = X(i,dim); %takes x1 value and puts into column
                        f(i,2) = i;
                    end
                    [val,idx]=sort(f(:,1));
                    f = [val idx];
                    
                    sp = [];
                    for i = 1:length(f)-1
                        
                        if (f(i,1) == f(i+1,1)) 
                            continue;    
                        end
                        sp = 0.5*(f(i,1) + f(i+1,1));
%                         left, right;
                        [w10, w20, lsum, rsum] = deal(0);
                        j = 1;
                        while (j <= i)
                            w10 = w10 + Y(f(j,2))*sv(f(j,2));
                            lsum = lsum + sv(f(j,2));
                            j = j+1;
                        end
                            
                        w10 = w10/lsum;
                        j = 1+i;
                        while (j<length(f)) %condition was as long as j is <= i
                            w20 = w20 + Y(f(j,2))*sv(f(j,2));
                            rsum = rsum + sv(f(j,2));
                            j = j+1;
                        end
                        w20 = w20/rsum;
                        
                        errl = 0; errr = 0;
                        j = 1;
                        while( j <= i)
                            errl = errl + (w10 - Y(f(j,2)))*(w10 - Y(f(j,2)))*sv(f(j,2));
                            j = j+1;
                        end
                        errl = errl/lsum;
                        
                        j = 1+i;
                        while(j <= length(f))
                            errr = errr + (w20 - Y(f(j,2)))*(w20 - Y(f(j,2)))*sv(f(j,2));
                            j = j + 1;
                        end
                        errr = errr/rsum;
                        
                        a = lsum/(lsum + rsum);
                        b = rsum/(lsum+rsum);
                        if (a*errl + b*errr < errBest || errBest == -1)
                            bestSplit = sp;
                            bestDim = dim;
                            errBest = a*errl + b*errr;
                            bestw10 = w10;
                            bestw20 = w20;
                        end
                    end
                end %best hardsplit ends ehre
                %(3) init params according to best hard split
                this.w = zeros(1:length(X(1,:)));
                for i = 1:length(this.w);
                    this.w(i) = 0.005+(-0.005-0.005).*rand(1,1);
                end
                this.w(bestDim) = -0.5;
                this.w0 = bestSplit*0.5;
%               // as described in the paper
%                 for i = 1:length(this.w);
%                     this.w(i) = 0.0;
%                 end
                
%                 this.w = this.w * 0;
%                 
%                 sigmoidSteepness = 10;
%                 
%                 this.w(bestDim) = sigmoidSteepness;
%                 this.w0 = -bestSplit*sigmoidSteepness;
                this.left.w0 = bestw10;
                this.right.w0 = bestw20;
        end %%end of hardinit
        
        function learnParams(this, X, Y, V, R, alpha, tree) 
        u = 0.1;
        eps = 0.00001;
        
        ix = [];
        
        dw = zeros(1:length(X(1,:))); %grads of w
        dwp = zeros(1:length(X(1,:))); %previous grads of 
        
        [dw10, dw20, dw0, dw10p, dw20p, dw0p] = deal(0.0);
        
        for i = 1:length(Y)
            ix(i) = i; 
        end
        MAXEPOCH = 25; 
        for e = 1:MAXEPOCH
            indx1 = randperm(length(ix));
            ix = ix(indx1);
%             ix = [13,19,14,15,12,4,18,3,2,0,9,16,10,11,6,17,7,1,5,8]+1;
            %for identical cpp implementation ^           
            for i = 1:length(X)
                j = ix(i);
                x = X(j,:); %should display both coordinates of corrseponding X
                r = Y(j); %target for this current X
                this.y = tree.evaluate(x); 
                d = this.y - r;
                
                t = alpha*d;
                m = this; 
                p = Node();
                while m.parent ~= 0
                    p = m.parent;
                    if m.isLeft
                        t = t * p.v;
                    else
                        t = t * (1-p.v);
                    end
                    m = m.parent;
                end
                
                dw = (-t * (this.left.y - this.right.y)*(this.v)*(1-this.v))*x;
                dw0 = -t * (this.left.y - this.right.y)*(this.v)*(1-this.v);
                dw10 = -t*(this.v);
                dw20 = -t*(1-this.v);
                
                %updating params:
                this.w = this.w + dw + u*dwp;
                this.w0 = this.w0 + dw0 + u*dw0p;
                this.left.w0 = this.left.w0 + dw10 + u*dw10p;
                this.right.w0 = this.right.w0 + dw20 + u*dw20p;    
                
                %update previous values
                dwp = dw;
                dw0p = dw0;
                dw10p = dw10;
                dw20p = dw20;
                
                alpha = alpha * 0.9999;
            end
        end
        end
        
        function this = splitNode(this, tree, X, Y, V, R)
        err = 0;
        type = tree.type;
        if (type == 'c')
            err = tree.errRate(V,R);
        else (type == 'r')
            err = tree.meanSqErr(V,R);
        end
        
        temp = this; %saving the current progress line 399(stores it at temp)
        save('temp_file.mat', 'temp');
        this.isLeaf = 0;
        this.w = zeros(1:length(X(1,:)));
        
        %left child
        this.left = Node();
        this.left.parent = this;
        this.left.isLeft = true;
        % right child
        this.right = Node();
        this.right.parent = this;
        this.right.isLeft = false;        
        
        bestErr = 1e10;
        MAXRETRY = 10;
        HARDINIT = true;%DUMMY VARIABLE CURRENTLY      
        %make MAXTRETY re-inits and choose the best
        for t = 1:MAXRETRY
            if(HARDINIT)
                this.hardInit(X,Y);
            else
                for i = 1:length(this.w)
                    this.w(i) = 0.005+(-0.005-0.005).*rand(1,1);
                end
                this.w0 = 0.005+(-0.005-0.005).*rand(1,1);
                this.left.w0 = 0.005+(-0.005-0.005).*rand(1,1);
                this.right.w0 = 0.005+(-0.005-0.005).*rand(1,1);
            end
            MINALPHA = 1;
            MAXALPHA = 10;
            alpha = MAXALPHA/(2^t);
            % alpha = 0; %mitigate learning!
            this.learnParams(X, Y, V, R, alpha, tree);
            
            if (type == 'r');
                newErr = tree.meanSqErr(V,R);
            else (type == 'c');
                newErr = tree.errRate(V,R);
            end
            if (newErr < bestErr)
                bestw = this.w;
                bestw0 = this.w0;
                bestw0l = this.left.w0;
                bestw0r = this.right.w0;
                bestErr = newErr;
            end
        end
        this.w = bestw;
        this.w0 = bestw0;
        this.left.w0 = bestw0l;
        this.right.w0 = bestw0r;
        PRETH = 1e-3;
        %continue splitting the children(457)
        if (bestErr + PRETH < err)
            this.left.splitNode(tree, X, Y, V, R);
            this.right.splitNode(tree, X, Y, V, R);
        else %stopping recursion?!
%             delete(this.left);
%             delete(this.right);
%             this.left = Node();
%             this.right = Node();
            load('temp_file.mat');
            this = temp;
%             this = temp.parent;
        end
        
        end
        
    end
    
end
