function savefile(Scalar,Pinf,Rho,I,II,...
    Td,Lg,Lc,xi,TIMESTEP,...
    B1RecData,B2RecData,B1Data,B2Data,WallData)

FLD.Pinf = Pinf;
FLD.Rho = Rho;

RecOption.Scalar = Scalar;
RecOption.ENDCYCLE = I;
RecOption.RecCYCLE = II;

Param.Td = Td;
Param.Lg = Lg;
Param.Lc = Lc;
Param.xi = xi;

filename = sprintf('');

save(filename,...
    'Param','FLD','RecOption','TIMESTEP',...
    'B1Data','B1RecData','B2RecData','B2Data','WallData')
end