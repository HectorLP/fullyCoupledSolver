//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////
Info<< "\nReading g for subsurface flow" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//info<< "\nReading g for subsurface flow" << endl;
uniformDimensionedVectorField gRef
(
    IOobject
    (
        "gRef",
        runTime.constant(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "\nReading hSubRef" << endl;
uniformDimensionedScalarField hSubRef
(
    IOobject
    (
        "hSubRef",
        runTime.constant(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("hSubRef", dimLength, 0)
);

dimensionedScalar ghRefSub
(
    mag(g.value()) > SMALL
    ? g & (cmptMag(g.value())/mag(g.value()))*hSubRef
    : dimensionedScalar("ghSub", g.dimensions()*dimLength, 0)
);
volScalarField ghSub("ghSub", (g & meshSubsurface.C()) - ghRefSub);
surfaceScalarField ghfSub("ghfSub", (g & meshSubsurface.Cf()) - ghRefSub);


word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Courant") CFLint = 1; 
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats and Courant available" 
            << exit(FatalError);
}

Info<< "Reading field pSub\n" << endl;
volScalarField pSub
(
    IOobject
    (
        "pSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
     meshSubsurface
);


Info<< "Reading field Sb \n" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Sbmin(transportProperties.lookupOrDefault("Sbmin",dimensionedScalar("Sbmin",dimless,0)));
dimensionedScalar Sbmax(transportProperties.lookupOrDefault("Sbmax",dimensionedScalar("Sbmax",dimless,1.0)));
dimensionedScalar D_sub(transportProperties.lookup("D_sub"));

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(meshSubsurface,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(meshSubsurface,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();	

volScalarField rhoSub
(
    IOobject
    (
        "rhoSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (1 - Sb) *rhoa + Sb*rhob
);
rhoSub.oldTime();

dimensionedScalar Mmu(mub/mua);//viscosity ratio

// Relative permeability model 
autoPtr<relativePermeabilityModel> krModel = relativePermeabilityModel::New("krModel",transportProperties,Sb);

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading Porous Media Properties \n" << endl;

// Porosity	
dimensionedScalar epsScalar(transportProperties.lookupOrDefault("eps",dimensionedScalar("",dimless,0.)));
   
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    meshSubsurface
);

// Intrinsic permeability       
Info<< "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

// interpolated permeability
surfaceScalarField Kf(fvc::interpolate(K,"K"));

Info<< "Reading porous transport properties" << endl;
dictionary porousTransport(transportProperties.subDict("porousTransport"));
word phaseName(porousTransport.lookupOrDefault<word>("phaseName","notDefined"));
dimensionedScalar Kd(porousTransport.lookupOrDefault("Kd",dimensionedScalar("Kd_default",dimensionSet(-1,3,0,0,0),0.)));
dimensionedScalar rs(porousTransport.lookupOrDefault("rs",dimensionedScalar("rs_default",dimensionSet(1,-3,0,0,0),0.)));
dimensionedScalar epsTotal(porousTransport.lookupOrDefault("epsTotal",dimensionedScalar("epsTotal_default",dimless,0.)));
dimensionedScalar lambda(porousTransport.lookupOrDefault("lambda",dimensionedScalar("lambda_default",dimensionSet(0,0,-1,0,0),0.)));

//Whether to solve the scalar C (non-reactive tracer)
Switch CSubTracerSwitch
(
    transportProperties.lookup("CSubTracerSwitch")
);

//Whether to solve reactive transport process in the subsurface
Switch ChemSubSwitch
(
    transportProperties.lookup("CheSubSwitch")
);

//Whether to include the pressure gradient of surface water to the flow calculation
Switch gradSurfPressureGradientSwitch
(
    transportProperties.lookup("gradSurfPressureGradientSwitch")
);

Info << nl << "porousTransport parameters" << nl << "{" << endl;
Info << "    " << "phaseName " << " : " << phaseName << endl;
Info << "    " << Kd.name() << " : " << Kd.value() << endl;
Info << "    " << rs.name() << " : " << rs.value() << endl;
Info << "    " << epsTotal.name() << " : " << epsTotal.value() << endl;
Info << "    " << lambda.name() << " : " << lambda.value() << endl;
Info << "} \n" << endl;

//- creation of dispersion model
autoPtr<dispersionModels> DeffModel= dispersionModels::New("DeffModel",transportProperties,meshSubsurface);

Info << "Computing retard coefficient R" << endl;
volScalarField R ("R", 1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps);
// volScalarField R
// (
//     IOobject
//     (

//     ),
//     1 + (1-epsTotal.value()) * rs.value() * Kd.value() / eps.value();
// );
/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading field USub\n" << endl;
volVectorField USub
(
    IOobject
    (
        "USub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

// #include "createPhi.H"
// surfaceScalarField phiP = phi;

Info<< "Reading/calculating face flux field phiSub\n" << endl;

surfaceScalarField phiSub
(
    IOobject
    (
        "phiSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(USub)
);
surfaceScalarField phiP = phiSub;
/////////////////////////////////////////////////////////////////////////////
////////////////////////// CAPILLARY PRESSURE ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));
autoPtr<capillarityModel> pcModel = capillarityModel::New("pcModel",transportProperties,Sb);

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        meshSubsurface,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    meshSubsurface,
    dimensionedScalar("",dimless/dimTime,0)
);

//- Looking for 2nd order time-scheme for source term
bool isBackwardScheme = false;
if (word(meshSubsurface.ddtScheme("source")) == "backward")
{
    isBackwardScheme = true;
}

Info<< "Reading field C\n" << endl;
volScalarField C_Sub
(
    IOobject
    (
        "C_Sub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

Info << "Reading solution component dictionary" << endl;
const dictionary& solutionComponentDict = transportProperties.subDict("solutionSpecies");
Info << "Reading solution component name\n" << endl;
speciesTable solutionSpecies
(
    solutionComponentDict.toc()
);
Info << "Reading each solution component's diffusion coefficient." << endl;
PtrList<dimensionedScalar> DY(solutionSpecies.size());

forAll(solutionComponentDict, i)
{
    const dictionary& subdict = solutionComponentDict.subDict(solutionSpecies[i]);
    DY.set
    (
        i, 
        new dimensionedScalar(subdict.lookup("D"))
    );
}


Info << "Create the fluid solution concentration field.\n" << endl;
PtrList<volScalarField> Y_Sub(solutionSpecies.size());

forAll(solutionSpecies, i)
{
    IOobject header
    (
        solutionSpecies[i],
        meshSubsurface.time().timeName(),
        meshSubsurface,
        IOobject::NO_READ
    );

    if (header.typeHeaderOk<volScalarField>(true))
    {
        Y_Sub.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    solutionSpecies[i],
                    meshSubsurface.time().timeName(),
                    meshSubsurface,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                meshSubsurface
            )
        );
    }
    else
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                meshSubsurface.time().timeName(),
                meshSubsurface,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            meshSubsurface
        );

        Y_Sub.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    solutionSpecies[i],
                    meshSubsurface.time().timeName(),
                    meshSubsurface,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Ydefault
            )
        );
    }
}

wordList surfaceSpeciesNamesSub;
if (transportProperties.found("surfaceSpecies"))
{
    Info << "Reading surface species dictionary" << endl;
    const dictionary& surfaceSpeciesDictSub = transportProperties.subDict("surfaceSpecies");
    Info << "Reading surface species name\n" << endl;
    surfaceSpeciesNamesSub = surfaceSpeciesDictSub.toc();
}

speciesTable surfaceSpeciesSub
(
    surfaceSpeciesNamesSub
);

Info << "Create surface composition vector\n" << endl;
PtrList<volScalarField> sYSub(surfaceSpeciesSub.size());

forAll(surfaceSpeciesSub, i)
{
    IOobject header
    (
        surfaceSpeciesSub[i],
        meshSubsurface.time().timeName(),
        meshSubsurface,
        IOobject::NO_READ
    );

    // check if field exists and can be read
    if (header.typeHeaderOk<volScalarField>(true))
    {
        sYSub.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    surfaceSpeciesSub[i],
                    meshSubsurface.time().timeName(),
                    meshSubsurface,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                meshSubsurface
            )
        );
    }
    else
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                meshSubsurface.time().timeName(),
                meshSubsurface,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            meshSubsurface
        );

        sYSub.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    surfaceSpeciesSub[i],
                    meshSubsurface.time().timeName(),
                    meshSubsurface,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Ydefault
            )
        );
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info << "Reading surface masters dictionary" << endl;
wordList surfaceMastersNamesSub;
if (transportProperties.found("surfaceMasters"))
{
    const dictionary& surfacesDictSub = transportProperties.subDict("surfaceMasters");
    Info << "Reading surface masters \n" << endl;
    surfaceMastersNamesSub = surfacesDictSub.toc();
}


speciesTable surfaceMastersSub
(
    surfaceMastersNamesSub
);

Info<< "Reading field Surf\n" << endl;
volScalarField SurfSub
(
    IOobject
    (
        "SurfSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface
);

Info << "Reading kinetic phases dictionary" << endl;
wordList kineticPhasesNamesSub;
if (transportProperties.found("kineticPhases"))
{
    Info << "Reading kinetic phases dictionary" << endl;
    const dictionary& kineticPhasesDictSub = transportProperties.subDict("kineticPhases");
    Info << "Reading kinetic phases names\n" << endl;
    kineticPhasesNamesSub = kineticPhasesDictSub.toc();
}

speciesTable kineticPhasesSub
(
    kineticPhasesNamesSub
);

PtrList<volScalarField> RSub(kineticPhasesSub.size());
forAll(kineticPhasesSub, i)
{
    RSub.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "RSub_" + kineticPhasesSub[i],
                meshSubsurface.time().timeName(),
                meshSubsurface,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            meshSubsurface,
            dimensionedScalar("RSub_"+kineticPhasesSub[i], dimensionSet(0,-2,-1,0,1,0,0), 0.0)
        )
    );
}

wordList selectedOutputNamesSub;
if (transportProperties.found("selectedOutput"))
{
    Info << "Reading selected output" << endl;
    const dictionary& selectedOutputDictSub = transportProperties.subDict("selectedOutput");
    Info << "Reading solution species name\n" << endl;
    selectedOutputNamesSub = selectedOutputDictSub.toc();
}


Info << "Create selected output vector\n" << endl;
PtrList<volScalarField> sOutSub(selectedOutputNamesSub.size());

forAll(selectedOutputNamesSub, i)
{
    sOutSub.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                selectedOutputNamesSub[i],
                meshSubsurface.time().timeName(),
                meshSubsurface,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            meshSubsurface
        )
    );
}

//surface potential
volScalarField psiSub
(
    IOobject
    (
        "psiSub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface,
    dimensionedScalar("psiSub", dimensionSet(1,2,-3,0,0,-1,0), 0.0)
);

//ionic strength
volScalarField ISub
(
    IOobject
    (
        "ISub",
        runTime.timeName(),
        meshSubsurface,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    meshSubsurface,
    dimensionedScalar("ISub", dimMoles/dimVolume, 0.0)
);

volScalarField alpha1Sub
(
    IOobject
    (
        "alpha1Sub",
        runTime.timeName(),
        meshSubsurface
    ),
    meshSubsurface,
    scalar(1.0)
);