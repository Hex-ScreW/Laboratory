%気泡クラス
classdef BubbleClass < handle
    properties
        % 計算用
        N = 60;         % 離散点数
        GaussNodeNum = 10; % わき出し間の分割数
        SCN = 0.2; % ソースを計算点から離す割合
        Node        % ルジャンドル求積 ノード
        Weight      % ルジャンドル求積 ウェイト
        Mip1 % 内挿関数
        Mip2
        
        % 気泡性質用
        MaxRadius = 0.0015; % 気泡の最大半径
        InitialRadius  % 気泡初期半径
        Radius
        Pc = 2300;   % Pa 気泡内の圧力
        RayleighTime % Rayleighの気泡崩壊時間
        MaxVolume % 最大体積
        
        % 座標
        X
        Y
        Xs
        Ys
        
        Q % 湧き出し量
        Phi % ポテンシャル
        Volume % 体積
        
        % 速度
        U
        V
    end
    methods
        %---------------コンストラクタ-------------------
        function obj = BubbleClass(N,Radius)
            
            if nargin > 0
                obj.N = N;
                obj.MaxRadius = Radius;
            end
            
            Initialize(obj);
        end
        
        function Initialize(obj)
            Initialdef(obj);
            CirclePlot(obj);
            PlotSource(obj);
            InitialPhi(obj);
            obj.Q = zeros(obj.N,1);
            obj.Volume = (4.*pi.*obj.InitialRadius.^3)./3;
            obj.U = zeros(1,obj.N);
            obj.V = zeros(1,obj.N);
        end
        
        function Initialize_g(obj)
            obj.InitialRadius = obj.MaxRadius;
            obj.Radius = obj.InitialRadius;
            obj.MaxVolume = (4.*pi.*obj.MaxRadius.^3)./3;
            obj.Volume = (4.*pi.*obj.Radius.^3)./3;
            obj.Phi = zeros(obj.N,1);
            obj.Pc = 101300;
            CirclePlot(obj);
        end
        
        %----------------------諸々の初期設定-------------------------------
        function Initialdef(obj)
            [obj.Node, obj.Weight] = GNAW(obj.GaussNodeNum); % ガウス求積のノードと重みを計算
            obj.Mip1 = 0.5.*(1-obj.Node); % ガウス求積の内挿関数
            obj.Mip2 = 0.5.*(1+obj.Node);
            obj.InitialRadius = obj.MaxRadius*0.05;
            obj.Radius = obj.InitialRadius;
            obj.RayleighTime = 1.8294*obj.MaxRadius*sqrt(998.203/(101300 - obj.Pc));
            obj.MaxVolume = 4.*pi.*(obj.MaxRadius.^3)./3;
        end
        
        %----------------------気泡表面プロット-------------------------------
        function CirclePlot(obj)
            Epsilon = 0.005/obj.N;
            theta = linspace(pi*0.5,-0.5*pi,obj.N);
            obj.Y = obj.Radius.*cos(theta);
            obj.X = obj.Radius.*sin(theta);
            
            %　北極の処理
            obj.Y(1,1) = obj.Radius.*cos(0.5*pi*(1-Epsilon));
            obj.X(1,1) = obj.Radius.*sin(0.5*pi*(1-Epsilon));
            %　南極の処理
            obj.Y(1,obj.N) = obj.Radius.*cos(-0.5*pi*(1-Epsilon));
            obj.X(1,obj.N) = obj.Radius.*sin(-0.5*pi*(1-Epsilon));
        end
        
        %------------------------ソース座標の計算-------------------------------
        function PlotSource(obj)
            %------表面要素の傾きの計算-----------
            DeltaY = circshift(obj.Y,-1) - circshift(obj.Y,1);
            DeltaX = circshift(obj.X,-1) - circshift(obj.X,1);
            
            TEMP = sqrt((obj.X' - obj.X).^2 + (obj.Y' - obj.Y).^2);
            TEMP(TEMP == 0) = 100;
            A = TEMP(end,1);
            TEMP = min(TEMP);
            B = TEMP(1);
            TEMP = TEMP.*obj.SCN;
            % 極が近づいたときに極間距離の1/3に修正
            if A == B
                A = A./3;
                TEMP(1) = A;
                TEMP(end) = A;
            end
            Distance = sqrt(DeltaX.^2+DeltaY.^2);
            dx = DeltaX./Distance.*TEMP;
            dy = DeltaY./Distance.*TEMP;
            
            obj.Xs = obj.X - dy;
            obj.Ys = obj.Y + dx;
            % 北極と南極のソースは別に処理する
            obj.Ys(1) = obj.Y(1);
            obj.Ys(obj.N) = obj.Y(obj.N);
            obj.Xs(1) = obj.X(1) - TEMP(1);
            obj.Xs(obj.N) = obj.X(obj.N) + TEMP(end);
        end
        
        %------------------------初期ポテンシャルの計算-------------------------------
        function InitialPhi(obj)
            obj.Phi = -obj.InitialRadius*sqrt((0.666667*((101300 - obj.Pc)/998.203)*...
                ((obj.MaxRadius/obj.InitialRadius)^3 - 1)));
            obj.Phi = ones(obj.N,1).*obj.Phi;
        end
        
        %--------------------要素間の最小距離を計算-------------------------
        function DISTANCE = MinimumElementDistance(obj)
            DISTANCE = sqrt((obj.X' - obj.X).^2 + (obj.Y' - obj.Y).^2);
            DISTANCE(DISTANCE == 0) = 1000;
            DISTANCE = min(DISTANCE(:));
        end
        
        %----------------------係数一括計算---------------------------
        function C = InteractionCoefficient(obj1, obj2, obj3)
            switch nargin
                case 1
                    C = Coefficient(obj1,obj1);
                case 2
                    C11 = Coefficient(obj1,obj1);
                    C12 = Coefficient(obj1,obj2);
                    C21 = Coefficient(obj2,obj1);
                    C22 = Coefficient(obj2,obj2);
                    C = [C11 C12; C21 C22];
                case 3
                    C11 = Coefficient(obj1,obj1);
                    C12 = Coefficient(obj1,obj2);
                    C13 = Coefficient(obj1,obj3);
                    C21 = Coefficient(obj2,obj1);
                    C22 = Coefficient(obj2,obj2);
                    C23 = Coefficient(obj2,obj3);
                    C = [C11 C12 C13; C21 C22 C23];
                otherwise
                    error('Please check the argment.')
            end
        end
        
        %------------------------係数の計算-----------------------------------
        % Cofficient('影響を受けるオブジェクト','影響を与えるオブジェクト')
        function C = Coefficient(obj1,obj2)
            % 計算点が1点のときのみ例外
            if obj2.N == 1
                R12 = (obj1.X' - obj2.Xs).^2 + (obj1.Y' - obj2.Ys).^2;
                R22 = (obj1.X' - obj2.Xs).^2 + (obj1.Y' + obj2.Ys).^2;
                k2 = 1 - R12./R22;
                Kk = ellipke(k2);
                C = -(obj2.Ys.*Kk)./(pi.*sqrt(R22));
            else
                % ガウス求積用座標
                XsZi = obj2.Mip1.*obj2.Xs + obj2.Mip2.*circshift(obj2.Xs,-1); XsZi(:,obj2.N) = [];
                YsZi = obj2.Mip1.*obj2.Ys + obj2.Mip2.*circshift(obj2.Ys,-1); YsZi(:,obj2.N) = [];
                % わき出し間の距離
                Li = sqrt(diff(obj2.Xs).^2 + diff(obj2.Ys).^2);
                for ZiN = 1:obj2.GaussNodeNum
                    % 1レイヤーを上下に2分割して計算
                    R12 = (obj1.X' - XsZi(ZiN,:)).^2 + (obj1.Y' - YsZi(ZiN,:)).^2;
                    R22 = (obj1.X' - XsZi(ZiN,:)).^2 + (obj1.Y' + YsZi(ZiN,:)).^2;
                    k2 = 1 - R12./R22;
                    Kk = ellipke(k2);
                    Greenfun = -(YsZi(ZiN,:).*Kk)./(pi.*sqrt(R22));
                    % S1の計算
                    S1 = Greenfun.*obj2.Mip1(ZiN).*obj2.Weight(ZiN);
                    S1 = S1.*Li.*0.5;
                    CS1 = [S1 zeros(obj1.N,1)];
                    % S2の計算
                    S2 = Greenfun.*obj2.Mip2(ZiN).*obj2.Weight(ZiN);
                    S2 = S2.*Li.*0.5;
                    CS2 = [zeros(obj1.N,1) S2];
                    % レイヤーの重ね合わせ
                    if ZiN == 1
                        C = CS1 + CS2;
                    else
                        C = C + CS1 + CS2;
                    end
                end
            end
        end
        
        %-------------------誘起速度一括計算----------------------
        function InducedVelocity(obj1,obj2,obj3)
            switch nargin
                case 1
                    [obj1.U,obj1.V] = CalVelocity(obj1,obj1);
                case 2
                    [U11,V11] = CalVelocity(obj1,obj1);
                    [U12,V12] = CalVelocity(obj2,obj1);
                    [U21,V21] = CalVelocity(obj1,obj2);
                    [U22,V22] = CalVelocity(obj2,obj2);
                    obj1.U = U11 + U12; obj1.V = V11 + V12;
                    obj2.U = U21 + U22; obj2.V = V21 + V22;
                case 3
                    [U11,V11] = CalVelocity(obj1,obj1);
                    [U12,V12] = CalVelocity(obj2,obj1);
                    [U13,V13] = CalVelocity(obj3,obj1);
                    [U21,V21] = CalVelocity(obj1,obj2);
                    [U22,V22] = CalVelocity(obj2,obj2);
                    [U23,V23] = CalVelocity(obj3,obj2);
                    [U31,V31] = CalVelocity(obj1,obj3);
                    [U32,V32] = CalVelocity(obj2,obj3);
                    [U33,V33] = CalVelocity(obj3,obj3);
                    obj1.U = U11 + U12 + U13; obj1.V = V11 + V12 + V13;
                    obj2.U = U21 + U22 + U23; obj2.V = V21 + V22 + V23;
                    obj3.U = U31 + U32 + U33; obj3.V = V31 + V32 + V33;
                otherwise
                    error('Please check the argment.')
            end
        end
        
        % 誘起速度の計算
        % CalVelocity('影響を与えるオブジェクト','影響を受けるオブジェクト')
        function [U, V] = CalVelocity(obj1,obj2)
            if obj1.N == 1
                R12 = (obj2.X - obj1.Xs).^2 + (obj2.Y - obj1.Ys).^2;
                R22 = (obj2.X - obj1.Xs).^2 + (obj2.Y + obj1.Ys).^2;
                k2 = 1 - R12./R22;
                [Kk,Ek] = ellipke(k2);
                %-----------------U速度成分--------------
                U = (obj2.X-obj1.Xs)./obj2.Y;
                U = U./sqrt(R22);
                U = U.*(k2./(1-k2));
                U = U.*Ek.*0.07957747;
                U = U.*obj1.Q;
                %----------------V速度成分----------------
                V = k2+(k2-2).*(obj1.Ys./obj2.Y);
                V = V./(1-k2);
                V = V.*Ek + 2.*(obj1.Ys./obj2.Y).*Kk;
                V = (V./sqrt(R22)).*0.07957747;
                V = V.*obj1.Q;
            else
                %　ガウス求積用座標
                XsZi = obj1.Mip1.*obj1.Xs + obj1.Mip2.*circshift(obj1.Xs,-1); XsZi(:,obj1.N) = [];
                YsZi = obj1.Mip1.*obj1.Ys + obj1.Mip2.*circshift(obj1.Ys,-1); YsZi(:,obj1.N) = [];
                %　わき出し間の距離
                Li = sqrt(diff(obj1.Xs).^2 + diff(obj1.Ys).^2);
                %　わき出しを分布化
                mZi = obj1.Mip1.*obj1.Q' + obj1.Mip2.*circshift(obj1.Q,-1)';
                mZi(:,obj1.N) = [];
                
                for ZiN = 1:obj1.GaussNodeNum
                    R12 = (obj2.X' - XsZi(ZiN,:)).^2 + (obj2.Y' - YsZi(ZiN,:)).^2;
                    R22 = (obj2.X' - XsZi(ZiN,:)).^2 + (obj2.Y' + YsZi(ZiN,:)).^2;
                    k2 = 1 - R12./R22;
                    [Kk,Ek] = ellipke(k2);
                    %-----------------U速度成分--------------
                    Uc = (obj2.X'-XsZi(ZiN,:))./obj2.Y';
                    Uc = Uc./sqrt(R22);
                    Uc = Uc.*(k2./(1-k2));
                    Uc = Uc.*Ek.*0.07957747;
                    Uc = Uc.*mZi(ZiN,:).*obj1.Weight(ZiN);
                    %----------------V速度成分----------------
                    Vc = k2+(k2-2).*(YsZi(ZiN,:)./obj2.Y');
                    Vc = Vc./(1-k2);
                    Vc = Vc.*Ek + 2.*(YsZi(ZiN,:)./obj2.Y').*Kk;
                    Vc = (Vc./sqrt(R22)).*0.07957747;
                    Vc = Vc.*mZi(ZiN,:).*obj1.Weight(ZiN);
                    
                    if ZiN == 1
                        U = Uc;
                        V = Vc;
                    else
                        U = U + Uc;
                        V = V + Vc;
                    end
                end
                U = U.*Li.*0.5;
                V = V.*Li.*0.5;
                U = sum(U,2)'; V = sum(V,2)';
            end
        end
        
        function temp = MaxMovement(obj,dt)
            dx = obj.X.*dt;
            dy = obj.Y.*dt;
            temp = sqrt(dx.^2 + dy.^2);
            temp = max(temp);
        end
        
        function UpdatePotential(obj,dt)
            obj.Phi = obj.Phi' + ((101300-obj.Pc)./998.203 + (obj.U.^2 + obj.V.^2).* 0.5).*dt;
            obj.Phi = obj.Phi';
        end
        
        %---------------------データを構造体に保存--------------------
        function s = saveobj(obj)
            s.N = obj.N;
            s.GaussNodeNum = obj.GaussNodeNum;
            s.SCN = obj.SCN;
            s.MaxRadius = obj.MaxRadius;
            s.InitialRadius = obj.InitialRadius;
            s.RayleighTime = obj.RayleighTime;
            s.Pc = obj.Pc;
        end
    end
    
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                newObj.N = s.N;
                newObj.GaussNodeNum = s.GaussNodeNum;
                newObj.SCN = s.SCN;
                newObj.MaxRadius = s.MaxRadius;
                newObj.InitialRadius = s.InitialRadius;
                newObj.RayleighTime = s.RayleighTime;
                newObj.Pc = s.Pc;
                obj = newObj;
            else
                obj = s;
            end
        end
    end
end