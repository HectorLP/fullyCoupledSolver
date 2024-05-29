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

Info<< "Reading field p\n" << endl;
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
