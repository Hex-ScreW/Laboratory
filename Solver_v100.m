close all
clear all
clc

maxNumCompThreads(2);
%maxNumCompThreads('automatic');

global Pinf Rho dt
Pinf = 101300;       % Pa�@�������ł̈���
Rho = 998.203;       % kg/m^3�@���̖��x�@20��,1atm
dt = 0.000000005;
ENDCYCLE = 200000;

RecOption.Scalar = 2;

II = 1;

TIMEtemp = 0.0;
TIME = zeros(1,ENDCYCLE/RecOption.Scalar);
TIMESTEP = zeros(1,ENDCYCLE/RecOption.Scalar);

%-----------------------���C���p�����[�^--------------------------------
Param.Td = 0.0; % �ʑ���
Param.Lg = 1.0; % �K�X�C�A�̕ǖʂ���̋���
Param.Lc = 1.5; % �L���r�e�[�V�����C�A�̕ǖʂ���̋���
Param.xi = 0.15; % �C�A���a��

%------------------�C�A�ݒ�------------------------
% �L���r�e�[�V�����C�A�̐ݒ�
Bc = BubbleClass(60,0.0015);

% �K�X�C�A�̐ݒ�
Bg = BubbleClass(1,Bc.MaxRadius.*Param.xi);
Bg.Initialize_g;
Bg.Y = Bg.Radius;

% �ǂ̐ݒ�
Wall = WallClass(24);

% �L�^�ݒ�
RecB1 = Record(Bc.N,ENDCYCLE/RecOption.Scalar);
RecB2 = Record(Bg.N,ENDCYCLE/RecOption.Scalar);
RecWall = Record(Wall.N,ENDCYCLE/RecOption.Scalar);

% �C�A�ԋ����̏C��
Bc.X = Bc.X - 2.*Param.Lc.*Bc.MaxRadius;
Bg.X = Bg.X - Bg.MaxRadius;

YCD(1) = 0.1*Bc.Y(1,2);
YCD(2) = 0.001*Bg.Y;

tic;
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
    if mod(I,RecOption.Scalar) == 0
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
    UpdatePotential(Bc);
    UpdatePotential(Bg);
    
    %---------�V�X�e���C���t�H���[�V����----------
    if rem(I,5000) == 0
        calculationtime = round(toc,1);
        fprintf('%g cycle %.1fsec\n',I,calculationtime);
    end
end
toc

RecOption.Pinf = Pinf;
RecOption.Rho = Rho;
RecOption.ENDCYCLE = I;
RecOption.RecCYCLE = II;

B1Data = Bc.saveobj;
B2Data = Bg.saveobj;
B1RecData = RecB1.saveobj(II,ENDCYCLE/RecOption.Scalar);
B2RecData = RecB2.saveobj(II,ENDCYCLE/RecOption.Scalar);
WallData = Wall.saveobj;
TIME(:,I:ENDCYCLE/RecOption.Scalar) = [];
TIMESTEP(:,I:ENDCYCLE/RecOption.Scalar) = [];
save('debug.mat','Param','B1Data','B1RecData','B2RecData','B2Data','WallData','RecOption','TIMESTEP')