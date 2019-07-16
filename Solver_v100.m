% �I�u�W�F�N�g�Ǘ� parfor�f�o�b�O

close all
clear all
clc

parfor PARCYCLE = 0:30
    try
        Pinf = 101300;       % Pa�@�������ł̈���
        Rho = 998.203;       % kg/m^3�@���̖��x�@20��,1atm
        dt = 0.000000005;
        ENDCYCLE = 200000;
        
        Scalar = 2;
        
        II = 1;
        
        TIMEtemp = 0.0;
        TIME = zeros(1,ENDCYCLE/Scalar);
        TIMESTEP = zeros(1,ENDCYCLE/Scalar);
        
        %-----------------------���C���p�����[�^--------------------------------
        Td = 0.0; % �ʑ���
        Lg = 1.0; % �K�X�C�A�̕ǖʂ���̋���
        Lc = 1.5; % �L���r�e�[�V�����C�A�̕ǖʂ���̋���
        xi = 0.1 + 0.01*PARCYCLE; % �C�A���a��
        
        %------------------�C�A�ݒ�------------------------
        % �L���r�e�[�V�����C�A�̐ݒ�
        Bc = BubbleClass(60,0.0015);
        
        % �K�X�C�A�̐ݒ�
        Bg = BubbleClass(1,Bc.MaxRadius.*xi);
        Bg.Initialize_g;
        Bg.Y = Bg.Radius;
        
        % �ǂ̐ݒ�
        Wall = WallClass(24,Bc.MaxRadius);
        
        % �L�^�ݒ�
        RecB1 = Record(Bc.N,ENDCYCLE/Scalar);
        RecB2 = Record(Bg.N,ENDCYCLE/Scalar);
        RecWall = Record(Wall.N,ENDCYCLE/Scalar);
        
        % �C�A�ԋ����̏C��
        Bc.X = Bc.X - 2.*Lc.*Bc.MaxRadius;
        Bg.X = Bg.X - Bg.MaxRadius;
        
        YCD = 0.1*Bc.Y(1,2);
        YCD(2) = 0.001*Bg.Y;
        
        for I = 1:ENDCYCLE
            %--------------------------�킫�o���_�̌���-----------------------------
            PlotSource(Bc);
            
            Bg.Xs = Bg.X;
            Bg.Ys = Bg.Y - Bg.Radius.*Bg.SCN;
            
            %-------------------------------�W��-----------------------------------
            CC = InteractionCoefficient(Bc,Bg,Wall);
            CV = InteractionCoefficient_from_Velocity(Wall,Bc,Bg);
            CMat = [CC ;CV];
            
            %---------------------�킫�o�������̌v�Z--------------------------------
            PhiMat = [Bc.Phi; Bg.Phi; Wall.Phi];
            QMat = CMat\PhiMat;
            Bc.Q = QMat(1:Bc.N);
            Bg.Q = QMat(Bc.N + 1:Bc.N + Bg.N);
            Wall.Q = QMat(Bc.N + Bg.N + 1:end);
            
            %--------------------------�U�N���x�̌v�Z-------------------------------
            InducedVelocity(Bc,Bg,Wall);
            
            %--------��ɖk�ɂ�V���폜-----------
            Bc.V(1,1) = 0.0; Bc.V(1,Bc.N) =  0.0;
            %     Bg.V(1,1) = 0.0; Bg.V(1,Bg.N) =  0.0;
            %     Bg.U(1,1) = 0.0;
            
            MaxSPEEDBc = max(sqrt(Bc.U.^2 + Bc.V.^2));
            MaxSPEEDBg = max(sqrt(Bg.U.^2 + Bg.V.^2));
            MaxSPEED = max([MaxSPEEDBc MaxSPEEDBg]);
            
            dt = (1/MaxSPEED)^2;
            if dt > 0.000000005
                dt = 0.000000005;
            elseif dt < 0.000000000001
                dt = 0.000000000001;
            end
            
            %     MaxMovement = max(SPEEDBg.*dt);
            %     if MaxMovement > Bg.Radius*0.5
            %         dt = MaxMovement./SPEEDBg;
            %     end
            
            TIMESTEP(I) = dt;
            TIMEtemp = TIMEtemp + dt;
            
            %---------------���W�̈ړ�----------------
            Bc.X = Bc.X + Bc.U.*dt;
            Bc.Y = Bc.Y + Bc.V.*dt;
            Bg.X = Bg.X + Bg.U.*dt;
            Bg.Y = Bg.Y + Bg.V.*dt;
            
            MinElemDist = [MinimumElementDistance(Bc) MinimumElementDistance(Bg)];
            MinElemDist = min(MinElemDist);
            
            %     if Bg.X(end) - Bg.X(1) > - 0.00005
            %         Bg.X(end) = Bg.X(1) - 0.00005;
            %         Bg.U(end) = 0;
            %     end
            
            %----------------�f�[�^�̋L�^--------------------
            if mod(I,Scalar) == 0
                RecordData(RecB1,Bc,II);
                RecordData(RecB2,Bg,II);
                RecordData(RecWall,Wall,II);
                TIME(II) = TIMEtemp;
                II = II + 1;
            end
            
            %-----------------�v���b�g�������щz���Ă邩����----------------
            Bc.Y = OverWriteY(Bc.Y,YCD(1));
            
            %----------------���͂̌v�Z------------------
            %Bg.Radius = abs(Bg.X(1) - Bg.X(end)).*0.5;
            Bg.Radius = Bg.Y;
            Bg.Pc = Pinf.*(Bg.MaxRadius./Bg.Radius).^3.9; % �̐ς���C�A�����͂��v�Z
            
            %-----------------���x�|�e���V�����̎��ԕω�---------------
            UpdatePotential(Bc,dt);
            UpdatePotential(Bg,dt);
            
            %---------�V�X�e���C���t�H���[�V����----------
            if rem(I,5000) == 0
                calculationtime = round(toc,1);
                fprintf('%g cycle %.1fsec\n',I,calculationtime);
            end
        end
    catch
        B1Data = Bc.saveobj;
        B2Data = Bg.saveobj;
        B1RecData = RecB1.saveobj(II,ENDCYCLE/Scalar);
        B2RecData = RecB2.saveobj(II,ENDCYCLE/Scalar);
        WallData = Wall.saveobj;
        TIME(:,I:ENDCYCLE/Scalar) = [];
        TIMESTEP(:,I:ENDCYCLE/Scalar) = [];
        
        savefile(Scalar,Pinf,Rho,I,II,...
            Td,Lg,Lc,xi,TIMESTEP,...
            B1RecData,B2RecData,B1Data,B2Data,WallData)
    end
end