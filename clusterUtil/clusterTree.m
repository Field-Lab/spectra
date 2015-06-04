classdef clusterTree < handle
    %CLUSTERTREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        root@clusterNode
        children@clusterTree;
    end
    
    methods
        function obj = clusterTree(low,high,pc)
            validateattributes(low,{'numeric'},{'integer','vector','nonempty','>',0},'','lower bound',1);
            validateattributes(high,{'numeric'},{'integer','vector','size',size(low),'>',0},'','upper bound',2);
            validateattributes(high-low,{'numeric'},{'vector','>=',0},'','cluster sizes');
%             validateattributes(order,{'numeric'},{'vector','integer','>',0},'','order',4);
%             validateattributes(projections,{'numeric'},{'2d','nrows',size(projections,1)},'','projections',3);
            
            obj.root = clusterNode(min(low),max(high),pc);
            
            argmaxes = intersect(find(high == max(high)),find(low == min(low)));
            if numel(argmaxes) > 0
                low(argmaxes(1)) = [];
                high(argmaxes(1)) = [];
            end
            for i = 1:numel(low)
                obj.insert(clusterTree(low(i),high(i),pc));
            end
        end
        
        function insert(obj,tree)
            % validateattributes(tree,{'clusterTree'},{'scalar'},'','tree',1);
            for k = 1:numel(obj.children)
                x = obj.children(k);
                if x.root.contains(tree.root)
                    x.insert(tree);
                    return
                end
                if tree.root.contains(x.root)
                    obj.children(k) = [];
                    tree.insert(x);
                    obj.insert(tree);
                    return
                end
            end
            
            for k = 1:numel(obj.children)
                x = obj.children(k);
                if x.root.overlap(tree.root)
                    tree.root = tree.root.merge(x.root);
                    obj.insert(tree)
                    return
                end
            end
            obj.children = [obj.children,tree];
        end
        
        function parentArray = assignPar(obj)
            parentArray = 0;
            
            for k = 1:numel(obj.children)
                x = obj.children(k).assignPar();
                parentArray = [parentArray,x+1];
            end
        end
        
        function reduce(obj) % Change to reduce (or) split single child into complementary
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
                nodeArray = [nodeArray,subNode];
                depth = [depth,subDepth+1];
            end
        end
        
        function [nodeArray,depth] = enumLeaves(obj)
            if numel(obj.children) == 0
                nodeArray = obj.root;
                depth = 1;
            else
                nodeArray = [];
                depth = [];
                for x = obj.children;
                    [subNode,subDepth] = x.enumLeaves();
                    nodeArray = [nodeArray,subNode];
                    depth = [depth,subDepth+1];
                end
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

