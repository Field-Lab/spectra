classdef clusterTree < handle
    %CLUSTERTREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        root@clusterNode
        children@clusterTree;
    end
    
    methods
        function obj = clusterTree(low,high,projections,order)
            validateattributes(low,{'numeric'},{'integer','vector','nonempty','>',0},'','lower bound',1);
            validateattributes(high,{'numeric'},{'integer','vector','size',size(low),'>',0},'','upper bound',2);
            validateattributes(high-low,{'numeric'},{'vector','>=',0},'','cluster sizes');
            validateattributes(order,{'numeric'},{'vector','integer','>',0},'','order',4);
            validateattributes(projections,{'numeric'},{'2d','nrows',size(projections,1)},'','projections',3);
            
            subproj = projections(order(min(low):max(high)),:);
            obj.root = clusterNode(min(low),max(high),mean(subproj,1),std(subproj,1));
            
            argmaxes = intersect(find(high == max(high)),find(low == min(low)));
            if numel(argmaxes) > 0
                low(argmaxes(1)) = [];
                high(argmaxes(1)) = [];
            end
            for i = 1:numel(low)
                obj.insert(clusterTree(low(i),high(i),projections,order));
            end
        end
        
        function insert(obj,tree)
           % validateattributes(tree,{'clusterTree'},{'scalar'},'','tree',1);
            found = false;
            for k = 1:numel(obj.children)
                x = obj.children(k);
                if x.root.contains(tree.root)
                    x.insert(tree);
                    found = true;
                    break
                end
                if tree.root.contains(x.root)
                    obj.children(k) = [];
                    tree.insert(x);
                    obj.insert(tree);
                    found = true;
                    break
                end
            end
            if ~found
                obj.children = [obj.children,tree];
            end
        end
        
        function parentArray = assignPar(obj)
            parentArray = 0;
            
            for k = 1:numel(obj.children)
                x = obj.children(k).assignPar();
                parentArray = [parentArray,x+1];
            end
        end
        
        function reduce(obj)
            if numel(obj.children) == 1
                obj.children = obj.children.children;
                obj.reduce();
            else
            for x = obj.children
                x.reduce();
            end
            end
        end
        
        function [nodeArray,depth] = enum(obj)
            nodeArray = obj.root;
            depth = 1;
            for x = obj.children;
                [subNode,subDepth] = x.enum();
                nodeArray = [subNode,nodeArray];
                depth = [subDepth+1,depth];
            end
        end
        
        function n = treeSize(obj)
            n = 1;
            for x = obj.children;
                n = n + x.treeSize;
            end
        end
    end
    
    
end

