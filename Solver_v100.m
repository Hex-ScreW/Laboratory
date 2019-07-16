close all
clear all
clc

maxNumCompThreads(2);
%maxNumCompThreads('automatic');

global Pinf Rho dt
Pinf = 101300;       % Pa　無限遠での圧力
Rho = 998.203;       % kg/m^3　水の密度　20℃,1atm
dt = 0.000000005;
ENDCYCLE = 200000;

RecOption.Scalar = 2;

II = 1;

TIMEtemp = 0.0;
TIME = zeros(1,ENDCYCLE/RecOption.Scalar);
TIMESTEP = zeros(1,ENDCYCLE/RecOption.Scalar);

%-----------------------メインパラメータ--------------------------------
Param.Td = 0.0; % 位相差
Param.Lg = 1.0; % ガス気泡の壁面からの距離
Param.Lc = 1.5; % キャビテーション気泡の壁面からの距離
Param.xi = 0.15; % 気泡半径比

%------------------気泡設定------------------------
% キャビテーション気泡の設定
Bc = BubbleClass(60,0.0015);

% ガス気泡の設定
Bg = BubbleClass(1,Bc.MaxRadius.*Param.xi);
Bg.Initialize_g;
Bg.Y = Bg.Radius;

% 壁の設定
Wall = WallClass(24);

% 記録設定
RecB1 = Record(Bc.N,ENDCYCLE/RecOption.Scalar);
RecB2 = Record(Bg.N,ENDCYCLE/RecOption.Scalar);
RecWall = Record(Wall.N,ENDCYCLE/RecOption.Scalar);

% 気泡間距離の修正
Bc.X = Bc.X - 2.*Param.Lc.*Bc.MaxRadius;
Bg.X = Bg.X - Bg.MaxRadius;

YCD(1) = 0.1*Bc.Y(1,2);
YCD(2) = 0.001*Bg.Y;

tic;
for I = 1:ENDCYCLE
    %--------------------------わき出し点の決定-----------------------------
    PlotSource(Bc);
    
    Bg.Xs = Bg.X;
    Bg.Ys = Bg.Y - Bg.Radius.*Bg.SCN;
    
    %-------------------------------係数-----------------------------------
    CC = InteractionCoefficient(Bc,Bg,Wall);
    CV = InteractionCoefficient_from_Velocity(Wall,Bc,Bg);
    CMat = [CC ;CV];
    
    %---------------------わき出し強さの計算--------------------------------
    PhiMat = [Bc.Phi; Bg.Phi; Wall.Phi];
    QMat = CMat\PhiMat;
    Bc.Q = QMat(1:Bc.N); 
    Bg.Q = QMat(Bc.N + 1:Bc.N + Bg.N);
    Wall.Q = QMat(Bc.N + Bg.N + 1:end);
    
    %--------------------------誘起速度の計算-------------------------------
    InducedVelocity(Bc,Bg,Wall);
    
    %--------南極北極のVを削除-----------
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
    
    %---------------座標の移動----------------
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
    
    %----------------データの記録--------------------
    if mod(I,RecOption.Scalar) == 0
        RecordData(RecB1,Bc,II);
        RecordData(RecB2,Bg,II);
        RecordData(RecWall,Wall,II);
        TIME(II) = TIMEtemp;
        II = II + 1;
    end
    
    %-----------------プロットが軸を飛び越えてるか判定----------------
    Bc.Y = OverWriteY(Bc.Y,YCD(1));
    
    %----------------圧力の計算------------------
    %Bg.Radius = abs(Bg.X(1) - Bg.X(end)).*0.5;
    Bg.Radius = Bg.Y;
    Bg.Pc = Pinf.*(Bg.MaxRadius./Bg.Radius).^3.9; % 体積から気泡内圧力を計算
    
    %-----------------速度ポテンシャルの時間変化---------------
    UpdatePotential(Bc);
    UpdatePotential(Bg);
    
    %---------システムインフォメーション----------
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