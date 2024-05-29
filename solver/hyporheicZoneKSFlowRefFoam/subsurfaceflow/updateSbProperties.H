//- relative permeability computation
Info << "Update the variables from the library." << endl;
krModel->correct();
Info << "Update kraf value." << endl;
kraf = fvc::interpolate(kra,"kra");
Info << "Update krbf value." << endl;
krbf = fvc::interpolate(krb,"krb");
Info << "Update dkrafdS value." << endl;
dkrafdS = fvc::interpolate(dkradS,"kra");
Info << "Update dkrbfdS value." << endl;
dkrbfdS = fvc::interpolate(dkrbdS,"krb");

//- mobility computation 
Info << "Update the value of Maf." << endl;
Maf = Kf*kraf/mua;
Info << "Update the value of Laf." << endl;
Laf = rhoa*Kf*kraf/mua;	
Info << "Update the value of Mbf." << endl;
Mbf = Kf*krbf/mub;
Info << "Update the value of Lbf." << endl;
Lbf = rhob*Kf*krbf/mub;
Info << "Update the value of Mf." << endl;
Mf = Maf+Mbf;
Info << "Update the value of Lf." << endl;
Lf = Laf+Lbf;
Info << "Update the value of Fbf." << endl;
Fbf = (krbf/mub) / ( (kraf/mua) + (krbf/mub) );

//- for source term in saturation equation
Info << "Update the value of Fb." << endl;
Fb = (krb/mub) / ( (kra/mua) + (krb/mub) );

//- capillarity computation
Info << "Update the capillarity value." << endl;
pcModel().correct();
// gradpc = fvc::interpolate(pcModel().dpcdS()*fvc::grad(Sb),"pc");

//- compute fluxes 
// phiPc = (Mbf & gradpc) & meshSubsurface.Sf();
Info << "Update the Pc flux." << endl;
phiPc = Mbf * fvc::interpolate(pcModel().dpcdS(),"pc")* fvc::snGrad(Sb) * meshSubsurface.magSf();
Info << "Update the flux with g." << endl;
phiG = (Lf * g) & meshSubsurface.Sf();
