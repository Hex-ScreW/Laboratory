%�ǃN���X
classdef WallClass < BubbleClass
    properties
        DeltaX
        DeltaY
    end
    methods
        %---------------�R���X�g���N�^-------------------
        function obj = WallClass(N,Radius)
            
            if nargin == 0
                TEMP = 24;
            else
                TEMP = N;
            end
            
            obj.N = TEMP;
            obj.MaxRadius = Radius;
            Initialize(obj);
            
        end
        
        function Initialize(obj)
            Initialdef(obj);
            InitialPlot(obj);
            PlotSource(obj);
            obj.Phi = zeros(obj.N,1);
            obj.Q = zeros(obj.N,1);
            obj.U = zeros(1,obj.N);
            obj.V = zeros(1,obj.N);
        end
        
        %----------------------�ǂ̏����v���b�g----------------------------
        function obj = InitialPlot(obj)
            %�@X,Y�����v���b�g�̌v�Z
            obj.Y = linspace(0,1,obj.N).^2;
            obj.Y = 6.*obj.MaxRadius.*obj.Y;
            obj.Y(1) = 1e-7;
            obj.X = zeros(1,obj.N);
            [obj.DeltaX, obj.DeltaY] = CalNormalVector(obj.X, obj.Y);
        end
        
        %---------------------�\�[�X���W�̌v�Z-----------------------------
        function obj = PlotSource(obj)
            obj.Xs = obj.X + obj.InitialRadius;
            obj.Ys = obj.Y;
        end
        
        %---------------------�ǖʗp�W���ꊇ�v�Z-----------------------------
        %�ǖʂƋC�A�̃C���^���N�V������O��ɂ��Ă���̂ŁA������2�`3
        function C = InteractionCoefficient_from_Velocity(obj1,obj2,obj3)
            if ~isa(obj1,'WallClass')
                error('Arg1 need to be WallClass. Please check the argment.')
            end
            switch nargin
                case 2
                    C1 = Coefficient_from_Velocity(obj1,obj2);
                    C2 = Coefficient_from_Velocity(obj1,obj1);
                    C = [C1 C2];
                case 3
                    C1 = Coefficient_from_Velocity(obj1,obj2);
                    C2 = Coefficient_from_Velocity(obj1,obj3);
                    C3 = Coefficient_from_Velocity(obj1,obj1);
                    C = [C1 C2 C3];
                otherwise
                    error('Please check the argment.')
            end
        end
        
        %-------------------------�ǖʗp�W���̌v�Z----------------------------
        %Coefficient_from_Velocity('�e�����󂯂�I�u�W�F�N�g�i�ǁj','�e����^����I�u�W�F�N�g')
        function C = Coefficient_from_Velocity(obj1,obj2)
            if ~isa(obj1,'WallClass')
                error('Arg1 need to be WallClass. Please check the argment.')
            end
            
            if obj2.N == 1
                R12 = (obj1.X' - obj2.Xs).^2 + (obj1.Y' - obj2.Ys).^2;
                R22 = (obj1.X' - obj2.Xs).^2 + (obj1.Y' + obj2.Ys).^2;
                k2 = 1 - R12./R22;
                [Kk,Ek] = ellipke(k2);
                %-----------------U���x����--------------
                Uc = (obj1.X'-obj2.Xs)./obj1.Y';
                Uc = Uc./sqrt(R22);
                Uc = Uc.*(k2./(1-k2));
                Uc = Uc.*Ek.*0.07957747;
                %----------------V���x����----------------
                Vc = k2+(k2-2).*(obj2.Ys./obj1.Y');
                Vc = Vc./(1-k2);
                Vc = Vc.*Ek + 2.*(obj2.Ys./obj1.Y').*Kk;
                Vc = (Vc./sqrt(R22)).*0.07957747;
                C = Uc + Vc;
            else
                %�@�K�E�X���ϗp���W
                XsZi = obj2.Mip1.*obj2.Xs + obj2.Mip2.*circshift(obj2.Xs,-1); XsZi(:,obj2.N) = [];
                YsZi = obj2.Mip1.*obj2.Ys + obj2.Mip2.*circshift(obj2.Ys,-1); YsZi(:,obj2.N) = [];
                %�@�킫�o���Ԃ̋���
                Li = sqrt(diff(obj2.Xs).^2 + diff(obj2.Ys).^2);
                for ZiN = 1:obj1.GaussNodeNum
                    %�@1���C���[���㉺��2�������Čv�Z
                    R12 = (obj1.X' - XsZi(ZiN,:)).^2 + (obj1.Y' - YsZi(ZiN,:)).^2;
                    R22 = (obj1.X' - XsZi(ZiN,:)).^2 + (obj1.Y' + YsZi(ZiN,:)).^2;
                    k2 = 1 - R12./R22;
                    [Kk,Ek] = ellipke(k2);
                    %-----------------U���x����--------------
                    Uc = (obj1.X'-XsZi(ZiN,:))./obj1.Y';
                    Uc = Uc./sqrt(R22);
                    Uc = Uc.*(k2./(1-k2));
                    Uc = Uc.*Ek.*0.07957747;
                    %----------------V���x����----------------
                    Vc = k2+(k2-2).*(YsZi(ZiN,:)./obj1.Y');
                    Vc = Vc./(1-k2);
                    Vc = Vc.*Ek + 2.*(YsZi(ZiN,:)./obj1.Y').*Kk;
                    Vc = (Vc./sqrt(R22)).*0.07957747;
                    %�@S1�̌v�Z
                    SU1 = Uc.*obj2.Mip1(ZiN).*obj2.Weight(ZiN);
                    SV1 = Vc.*obj2.Mip1(ZiN).*obj2.Weight(ZiN);
                    SU1 = SU1.*Li.*0.5;
                    SV1 = SV1.*Li.*0.5;
                    CSU1 = [SU1 zeros(obj1.N,1)];
                    CSV1 = [SV1 zeros(obj1.N,1)];
                    %�@S2�̌v�Z
                    SU2 = Uc.*obj2.Mip2(ZiN).*obj1.Weight(ZiN);
                    SV2 = Vc.*obj2.Mip2(ZiN).*obj1.Weight(ZiN);
                    SU2 = SU2.*Li.*0.5;
                    SV2 = SV2.*Li.*0.5;
                    CSU2 = [zeros(obj1.N,1) SU2];
                    CSV2 = [zeros(obj1.N,1) SV2];
                    %�@���C���[�̏d�ˍ��킹
                    if ZiN == 1
                        CCU = CSU1 + CSU2;
                        CCV = CSV1 + CSV2;
                    else
                        CCU = CCU + CSU1 + CSU2;
                        CCV = CCV + CSV1 + CSV2;
                    end
                end
                C = CCU.*obj1.DeltaY' + CCV.*obj1.DeltaX';
            end
        end
        
        %---------------------�f�[�^���\���̂ɕۑ�--------------------
        function s = saveobj(obj)
            s.N = obj.N;
            s.GaussNodeNum = obj.GaussNodeNum;
            s.SCN = obj.SCN;
            s.MaxRadius = obj.MaxRadius;
            s.InitialRadius = obj.InitialRadius;
            s.Pc = obj.Pc;
            s.X = obj.X;
            s.Y = obj.Y;
            s.Xs = obj.Xs;
            s.Ys = obj.Ys;
        end
    end
end