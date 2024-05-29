// --- Coats (2003) IMPES STABILITY: CFL LIMIT
{
    // - inertial part
    volScalarField dfw (
        "dfw",
        (dkrbdS*kra - dkradS*krb) /(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2))
    );
    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    // - gravity part
    dfw -= K*(rhoa-rhob)*fvc::surfaceSum(mag(meshSubsurface.Sf() & g))/fvc::surfaceSum(mag(phi)+smallRate)* (Foam::pow(kra,2)*dkrbdS/mua + Foam::pow(krb,2)*dkradS/mub)/(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2));;

    scalarField CFLCoats((runTime.deltaT()/eps)*dfw*fvc::surfaceSum(mag(phi)));

    // - capillarity part of CFL
    if(activateCapillarity)
    {
        CFLCoats += (runTime.deltaT()/eps)*2*mag(pcModel->dpcdS())*fvc::surfaceSum(Kf*meshSubsurface.magSf()/mag(meshSubsurface.delta()))*(kra*krb/(mub*kra+mua*krb));
    }

    CFLCoats /= meshSubsurface.V();
    CFLUse = gMax(cmptMax(CFLCoats));
    maxDeltaTFact = maxCo/(CFLUse + SMALL);

    Info<< "Coats CFL Number max: " << CFLUse << endl;
}
