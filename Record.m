classdef Record < handle
    properties
        X
        Y
        Xs
        Ys
        Phi
        Q
        U
        V
        Speed
    end
    methods
        function obj = Record(N, ENDCYCLE)
            if nargin == 2
                obj.X = zeros(N,ENDCYCLE);
                obj.Y = zeros(N,ENDCYCLE);
                obj.Xs = zeros(N,ENDCYCLE);
                obj.Ys = zeros(N,ENDCYCLE);
                obj.Phi = zeros(N,ENDCYCLE);
                obj.Q = zeros(N,ENDCYCLE);
                obj.U = zeros(N,ENDCYCLE);
                obj.V = zeros(N,ENDCYCLE);
            end
        end
        
        function RecordData(Recobj,obj,I)
            Recobj.X(:,I) = obj.X';
            Recobj.Y(:,I) = obj.Y';
            Recobj.Xs(:,I) = obj.Xs';
            Recobj.Ys(:,I) = obj.Ys';
            Recobj.Phi(:,I) = obj.Phi';
            Recobj.Q(:,I) = obj.Q';
            Recobj.U(:,I) = obj.U';
            Recobj.V(:,I) = obj.V';
        end
        
        function s = saveobj(obj, I, ENDCYCLE)
            obj.X(:,I:ENDCYCLE) = [];
            obj.Y(:,I:ENDCYCLE) = [];
            obj.Xs(:,I:ENDCYCLE) = [];
            obj.Ys(:,I:ENDCYCLE) = [];
            obj.Phi(:,I:ENDCYCLE) = [];
            obj.Q(:,I:ENDCYCLE) = [];
            obj.U(:,I:ENDCYCLE) = [];
            obj.V(:,I:ENDCYCLE) = [];
            obj.Speed = sqrt(obj.U.^2 + obj.V.^2);
            s.X = obj.X;
            s.Y = obj.Y;
            s.Xs = obj.Xs;
            s.Ys = obj.Ys;
            s.Phi = obj.Phi;
            s.Q = obj.Q;
            s.U = obj.U;
            s.V = obj.V;
            s.Speed = obj.Speed;
        end
    end
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                newObj.X = s.X;
                newObj.Y = s.Y;
                newObj.Phi = s.Phi;
                newObj.Q = s.Q;
                newObj.U = s.U;
                newObj.V = s.V;
                newObj.Speed = Speed.V;
                obj = newObj;
            else
                obj = s;
            end
        end
    end
end