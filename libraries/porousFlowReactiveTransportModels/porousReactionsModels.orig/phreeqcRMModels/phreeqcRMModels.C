/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "phreeqcRMModels.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "subCycleTime.H"
#include "subCycle.H"

namespace Foam
{
    namespace porousReactionsModels
    {
        defineTypeNameAndDebug(phreeqcRMModels,0);

        addToRunTimeSelectionTable
        (
            basicPorousReactionsModels,
            phreeqcRMModels,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousReactionsModels::phreeqcRMModels::phreeqcRMModels
(
    const fvMesh& mesh,
    const dictionary& dict 
)
:
    basicPorousReactionsModels(mesh, dict),
    phreeqcDict_(dict.subDict(typeName)),
    transportPropertiesDict_(dict),
    activateUpdatePorosity_(phreeqcDict_.lookup("activateUpdatePorosity")),
    setComponentH2O_
    (
            phreeqcDict_.lookupOrDefault<Switch>("setComponentH2O", false)
    ),
    numThreads_(phreeqcDict_.lookupOrDefault("numThreads", 1)),
    nthread_(1),
    nxyz_(mesh_.cells().size()),
    strangSteps_(readScalar(phreeqcDict_.lookup("StrangSteps"))),
    phreeqcInputFile_(phreeqcDict_.lookup("PhreeqcInputFile")),
    phreeqcDataBase_(phreeqcDict_.lookup("PhreeqcDataBase")),
    mineralSubDict_(mineralList_.size()),
    activatePhaseEquilibrium_(mineralList_.size()),
    Vm_(mineralList_.size()),
    linkDomainZonesToPhreeqcSol_
    (
        IOobject
        (
            "linkDomainZonesToPhreeqcSol",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        1
    ),
    cvODE (phreeqcDict_.lookupOrDefault("cvODE", false)),
    cvODETol_ (phreeqcDict_.lookupOrDefault("cvODETol", 1e-8)),
    use_SNIA (phreeqcDict_.lookupOrDefault("use_SNIA", false)),
    UName_(phreeqcDict_.lookupOrDefault<word>("U","U")),
    U_(mesh.lookupObject<volVectorField>(UName_)),
    pH_ 
    (
        IOobject
        (
            "pH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pH", dimless, 0.0),
        "zeroGradient"
    ),
    phreeqc_(mesh_.cells().size(), nthread_)
    {
        Info << "-------------------------------------------------------------" << endl
		   << "-------- Initialization of geochemistry with PhreeqC --------"
		   << nl << endl;
        phreeqc_.LoadDatabase(mesh_.time().constant()/phreeqcDataBase_);

        Info << "Set properties...";
    
        status = phreeqc_.SetErrorHandlerMode(1);
        status = phreeqc_.SetComponentH2O(setComponentH2O_);
        status = phreeqc_.SetRebalanceFraction(0.5);
        status = phreeqc_.SetRebalanceByCell(true);
        phreeqc_.UseSolutionDensityVolume(false);
        phreeqc_.SetPartitionUZSolids(false);
        Info << "OK" << nl << endl;

        Info << "Set concentration units... ";

        status = phreeqc_.SetUnitsSolution(2);  //1.mg/L; 2. mol/L; 3. kg/kgs
        status = phreeqc_.SetUnitsPPassemblage(1);  //0. mol/L cell; 1. mol/L water; 2. mol/L roc
        status = phreeqc_.SetUnitsExchange(1);  //0. mol/L cell; 1. mol/L water; 2. mol/L rock
        status = phreeqc_.SetUnitsSurface(1);   //0. mol/L cell; 1. mol/L water; 2. mol/L rock
        status = phreeqc_.SetUnitsGasPhase(1);  //0. mol/L cell; 1. mol/L water; 2. mol/L rock
        status = phreeqc_.SetUnitsSSassemblage(1); //0. mol/L cell; 1. mol/L water; 2. mol/L rock
        status = phreeqc_.SetUnitsKinetics(1);  //0. mol/L cell; 1. mol/L water; 2. mol/L rock

        Info << "OK" << nl << endl;

        Info << "Set volume, initial porosity and saturation...";
        std::vector<double> rv, por, sat;
        rv.resize(nxyz_,1.0);
        por.resize(nxyz_, 0);
        sat.resize(nxyz_, 1.0);

        forAll(eps_, cellI)
        {
            por[cellI] = eps_[cellI];

        }
        status = phreeqc_.SetRepresentativeVolume(rv);
        status = phreeqc_.SetPorosity(por);
        status = phreeqc_.SetSaturation(sat);

        Info << "OK" << nl << endl;

        bool workers = true;    //worker instances do the reaction calculation for transport 
        bool initial_phreeqc = true;    //InitialPhreeqc instance accumulate initial and boundary conditions
        bool utility = true;    //Utility instance is available for processing

        status = phreeqc_.SetSpeciesSaveOn(true);

        status = phreeqc_.RunFile
        (
            workers, initial_phreeqc, utility, mesh_.time().constant()/phreeqcInputFile_
        );
        //Determine number of the reaction module information
        int ncomps = phreeqc_.FindComponents();

        Info << 'Number of components = ' << ncomps << nl << endl;

        {
            Info << "Database:                                         "
                << phreeqc_.GetDatabaseFileName().c_str() << "\n";
            
            Info << "Number of threads:                                "
                << phreeqc_.GetThreadCount() << "\n";
            
            Info << "Number of MPI processes:                          "
                << phreeqc_.GetMpiTasks() << "\n";
            
            Info << "MPI task number:                                  "
                << phreeqc_.GetMpiMyself() << "\n";
            
            Info << "File prefix:                                      "
                << phreeqc_.GetFilePrefix() << "\n";
            
            Info << "Number of grid cells in the user's model:         "
                << phreeqc_.GetGridCellCount() << "\n";

            Info << "Number of chemistry cells in the reaction module: "
                << phreeqc_.GetChemistryCellCount() << "\n";

            Info << "Number of component for transport:                "
                << phreeqc_.GetComponentCount() << "\n";

            Info << "Partioning of UZ solids:                          "
                << phreeqc_.GetPartitionUZSolids() << "\n";

            Info << "Error handler mode:                               "
                << phreeqc_.GetErrorHandlerMode() << "\n" << endl;
        } 
        activateSelectedOutput();

        initializeMineralDistribution();

        initializeFluidComposition();

        if (status < 0) phreeqc_.DecodeError(status);

        Info << "\n--------------- End initialization of PhreeqC ---------------" << endl
                << "-------------------------------------------------------------"
                << nl << endl;

}

std::string Foam::porousReactionsModels::phreeqcRMModels::generateEqPhasesInputString()
{
    std::string input;

    forAll(mesh_.C(), i)
    {
        input += 
                " EQUILBRIUM_PHASES " + std::to_string(i) + "; \n";
        
        forAll(mineralList_, s)
        {
            if(activatePhaseEquilibrium_[s] == true)
            {
                double mineralMoleI = Ys_[s][i]/Vm_[s].value()*1e-3/(eps_[i]+VSMALL);
                std::ostringstream strs;
                strs << mineralMoleI;
                input += mineralList_[s] + "0.0 " + strs.str() + ";\n"; 
            }
        }
        
        input += "SAVE SOLUTION " + std::to_string(i) + " ; ";
        
        input += " END; \n";
        
    }

    return input;
}

std::string Foam::porousReactionsModels::phreeqcRMModels::generateKineticsInputString()
{
    std::string input;
    
    forAll(mesh_.C(), i)
    {
        input += " KINETICS " + std::to_string(i) + " ;\n ";

        forAll(mineralList_, s)
        {
            if(activatePhaseEquilibrium_[s] == false)
            {
                double mineralMoleI = Ys_[s][i]/Vm_[s].value() * 1e-3/(eps_[i]+VSMALL);
                std::ostringstream strs_Mi;
                strs_Mi << mineralMoleI; 

                const volScalarField Ae_ ("Ae", porousMedia_[s].surfaceArea());

                double AeMi = Ae_[i]/(Ys_[s][i]+SMALL);

                std::ostringstream strs_AeMi;
                strs_AeMi << AeMi;
                if (cvODE == false)
                {
                    Info << "Using RK to solve chemistry" << nl;
                    input += mineralList_[s] + ";\n"
                            +" -m " + strs_Mi.str() + ";\n"
                            +" -m0 " + strs_Mi.str() + ";\n"
                            +" -parms " + strs_AeMi.str() + " 0.0 " + ";\n";
                }
                else
                {
                    Info << "Using cv ode for chemistry " << nl;
                    std::ostringstream strs_cvODETol;
                    strs_cvODETol << cvODETol_;
                    input += mineralList_[s] + ";\n"
                        +" -m " + strs_Mi.str() + ";\n"
                        +" -m0 " + strs_Mi.str() + ";\n"
                        +" -parms " + strs_AeMi.str() + " 0.0 " + ";\n"
                        +" -cvode " + "true" + ";\n"
                        +" -cvode_steps " + "500" + ";\n"
                        +" -tol " + strs_cvODETol.str() + "; \n";
                }
            }
            
        } 
        input += " END; \n";
    }
    return input;
}

void Foam::porousReactionsModels::phreeqcRMModels::updateKineticsParameters()
{
    std::string input;

    forAll(mesh_.C(),i)
    {
        input += " KINETICS_MODIFY " + std::to_string(i) + " ; \n";
        forAll(mineralList_,s)
        {
            if(activatePhaseEquilibrium_[s] == false)
            {
                const volScalarField Ae_ ("Ae",porousMedia_[s].surfaceArea());
                //attention au 100 pour Calcite mais pas pour le reste
                double AeMi;
                if (eps_[i]<=0.001)
          		{
          			AeMi = SMALL;
          		}
          		else
          		{
          			AeMi = Ae_[i];
          		}
          		/*Saideep: to curtail add. precip. at low porosity*/
          		std::ostringstream strs_AeMi;
          		strs_AeMi << AeMi;
              input +=
              + "-component " + mineralList_[s]	+ ";\n"
              +" -d_params "  + strs_AeMi.str() + " 0.0" + ";\n"
              ;
            }
            input += " END; \n";
        }
  }
  status = phreeqc_.RunString(true, true, true, input.c_str());
}

void Foam::porousReactionsModels::phreeqcRMModels::initializeMineralDistribution()
{
    forAll(mineralList_, s)
    {
        word currentMineral = mineralList_[s];
        Info << " Doing stuff for mineral: " << currentMineral << endl;
        mineralSubDict_.set
        (
            s,
            new dictionary
            (
                transportPropertiesDict_.subDict(currentMineral+"Properties")
            )
        );
        activatePhaseEquilibrium_.set
        (
            s,
            new Switch
            (
                mineralSubDict_[s].lookup("activatePhaseEquilibrium")
            )
        );
        Vm_.set
        (
            s,
            new dimensionedScalar
            (
                mineralSubDict_[s].lookup("Vm")
            )
        );
    }

    forAll(mineralList_, s)
    {
        if(activatePhaseEquilibrium_[s] == true)
        {
            Info << mineralList_[s] << " is PHASE_EQ "
                << endl;
        }
        if(activatePhaseEquilibrium_[s] == false)
        {
            Info << mineralList_[s]
                << " is KINETICS "
                << endl;
        }
    }

    std::string input;
    input += generateEqPhasesInputString();
    input += generateKineticsInputString();

    status = phreeqc_.RunString(true, true, true, input.c_str());

    updatePorosityPhreeqc();
}

void Foam::porousReactionsModels::phreeqcRMModels::updatePorosityPhreeqc()
{
    if(activateUpdatePorosity_)
    {
        eps_ = 0.0 * eps_;
        forAll(mineralList_, s)
        {
            eps_ += Ys_[s];
        }
        eps_ = 1. - eps_ - inertMineral_;
        eps_.max(1e-3);
        
        std::vector<double> por;
        por.resize(nxyz_, 0);

        forAll(eps_, cellI)
        {
            por[cellI] = eps_[cellI];
        }

        status = phreeqc_.SetPorosity(por);
    }
}

void Foam::porousReactionsModels::phreeqcRMModels::updatePorosity()
{
    updatePorosityPhreeqc();
}

void Foam::porousReactionsModels::phreeqcRMModels::initializeFluidComposition()
{
    const std::vector<std::string> & components = phreeqc_.GetComponents();

    //Determine species information
    const std::vector<std::string> &species = phreeqc_.GetSpeciesNames();
    const std::vector<double> &species_z = phreeqc_.GetSpeciesZ();

    const std::vector<double> &species_d = phreeqc_.GetSpeciesD25();
    bool species_on = phreeqc_.GetSpeciesSaveOn();
    if (species_on == true) 
    {
        Info << "Species Save on " << nl << endl;
    }
    int nspecies = phreeqc_.GetSpeciesCount();

    Info << "nspecies = " << nspecies << endl;

    for (int i = 0; i < nspecies; ++i)
    {
        Info << species[i] << "\n";
        Info << "     Charge: " << species_z[i] << endl;
        Info << "     Dw:     " << species_d[i] << endl;
    }
    Info << endl;
	std::vector<int> ic1;
	ic1.resize(nxyz_*7, -1);
    for (int i = 0; i < nxyz_; ++i)
    {
        ic1[0*nxyz_+i] = linkDomainZonesToPhreeqcSol_[i];
        ic1[1*nxyz_+i] = i;
        ic1[2*nxyz_+i] = -1;
        ic1[3*nxyz_+i] = -1;
        ic1[4*nxyz_+i] = -1;
        ic1[5*nxyz_+i] = -1;
        ic1[6*nxyz_+i] = i;
    }
    status = phreeqc_.InitialPhreeqc2Module(ic1);

    //Initial equilibration of cells
    double time = 0.0;
    double time_step = 0.0;
    std::vector<double> c;
    c.resize(nxyz_*components.size());
    status = phreeqc_.SetTime(time);
    status = phreeqc_.SetTimeStep(time_step);
    status = phreeqc_.RunCells();
    status = phreeqc_.GetConcentrations(c);

    Info << "Get boundary conditions..." << endl;

    std::vector<double> bc_conc;
    std::vector<int> bc1;
    int nbound = 1;
    bc1.resize(nbound, 0);
    phreeqc_.InitialPhreeqc2Concentrations(bc_conc, bc1);
    Info << "OK" << nl << endl;

    //Transfer to the Unit of OpenFOAM 
    Info << "Set concentrations fields and transport properties with OpenFOAM....." << endl;

    int ncomps = phreeqc_.FindComponents();
    ncomps = phreeqc_.FindComponents();
    wordList componentNames(ncomps);
    forAll(componentNames, i)
    {
        componentNames[i] = components[i];
    }

    Y_.resize(componentNames.size());

    //-----Define boundary condition type
    label patchInlet = mesh_.boundaryMesh().findPatchID("inlet");
    if(patchInlet < 0)
    {
        Info << "ERROR: \"inlet\" is not defined as a boundary condition "
                <<  nl << endl;
    }
    wordList YiBoundaryType (mesh_.boundaryMesh().size(),"zeroGradient");
    forAll(YiBoundaryType,typeID)
    {
        if(mesh_.boundaryMesh().types()[typeID]=="empty")
        {
            YiBoundaryType[typeID]="empty";
        }
    }
    YiBoundaryType[patchInlet]="fixedValue";

    // ----- Create and initialize mixture
    forAll(componentNames,s)
    {
        word currentSpecies = componentNames[s];
        Info << " Doing stuff for components: " << currentSpecies << endl;

        Y_.set
        (
            s,
            new volScalarField
            (
                IOobject
                (
                    "Y."+componentNames[s],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(currentSpecies,dimless,0.0),
                YiBoundaryType
            )
        );

        //BC
        Y_[s].boundaryFieldRef()[patchInlet] == bc_conc[s];

        forAll(Y_[s],cellI)
        {
            Y_[s][cellI]=c[nxyz_*s+cellI];
        }

        Y_[s].write();

    }
    updatepH();
    pH_.write();

    Info << "OK"<< nl <<endl;
}

void Foam::porousReactionsModels::phreeqcRMModels::updateFluidComposition()
{

    if(use_SNIA)
    {
        Info << "Using SNIA (seq. non-iterative) for transport and chemistry coupling." << nl;

        const volTensorField &Deff = effectiveDispersionTensor();
        word divPhiYiScheme = "div(phi,Yi)";

        std::vector<double> c;
        c.resize(Y_.size()*nxyz_);

        forAll(Y_,i)
        {
             volScalarField& Yi = Y_[i];
             fvScalarMatrix YiEqn
             (
                 fvm::ddt(eps_,Yi) + fvm::div(phi_,Yi,divPhiYiScheme)
                 -fvm::laplacian(eps_*Deff,Yi,"laplacian(Di,Yi)")
             );
             YiEqn.solve();

             forAll(Y_[i],cellI)
             {
                 c[i*nxyz_+cellI]=Y_[i][cellI];
             } //copy flow information to phreeqc concentration
        } //solve transport for the species    i
        updateKineticsParameters();
        Info << "run Phreeqc ... " << nl;
        status = phreeqc_.SetConcentrations(c);
        status = phreeqc_.SetTimeStep(mesh_.time().deltaT().value());
        status = phreeqc_.SetTime(mesh_.time().value());
        Info << "I am going to runcell" << nl;
        status = phreeqc_.RunCells();
        Info << "I finished runcell" << nl;
        status = phreeqc_.GetConcentrations(c);

        forAll(Y_,s)
        {
            forAll(Y_[s],cellI)
            {
                Y_[s][cellI]=c[nxyz_*s+cellI];
            }
        }
        updatepH();
    }
    else
    {

        int strangCounter_ = 0;
  //    const volScalarField &Deff = effectiveDispersion();
        const volTensorField &Deff = effectiveDispersionTensor();


        word divPhiYiScheme = "div(phi,Yi)";

        std::vector<double> c;
        c.resize(Y_.size()*nxyz_);
        const Time& runTime = mesh_.time();

        for
        (
            subCycleTime subCycle(const_cast<Time&>(runTime),strangSteps_); !(++subCycle).end();
        )
        {
            strangCounter_++;
            Info << "Num sub cycles = " << subCycle.nSubCycles() << endl; //report numSubCycles
            Info << "Delta t = " << runTime.deltaT() << endl;             //subCycle time = dt/nSubCycles

            forAll (Y_, i)
            {
                volScalarField& Yi = Y_[i];
                fvScalarMatrix YiEqn
                (
                    fvm::ddt(eps_,Yi) + fvm::div(phi_,Yi,divPhiYiScheme)
                    -fvm::laplacian(eps_*Deff,Yi,"laplacian(Di,Yi)")
                );
                YiEqn.solve();

                forAll(Y_[i],cellI)
                {
                    c[i*nxyz_+cellI]=Y_[i][cellI];
                } //copy flow information to phreeqc concentration
            } //solve transport for the species

            for (; (strangCounter_==strangSteps_/2 /*|| strangCounter_==7*/); strangCounter_++)
            {
                Info<<"run phreeqc ..." << strangCounter_ << endl;
                status = phreeqc_.SetConcentrations(c);         // Transported concentrations

                status = phreeqc_.SetTimeStep(runTime.deltaT().value()*strangSteps_);    // Time step for kinetic reactions
                Info << "Time stamp = " << mesh_.time().timeName() << endl;
                status = phreeqc_.RunCells();
                // Transfer data from PhreeqcRM for transport
                status = phreeqc_.GetConcentrations(c);
                forAll (Y_,s)
                {
                    forAll(Y_[s],cellI)
                    {
                        Y_[s][cellI]=c[nxyz_*s+cellI];
                    }
                }
                phreeqc_.CloseFiles();
                //if (strangCounter_==3){
                //updateMineralDistribution();
                //updatePorosity();
                updateKineticsParameters();//}
                updatepH();
            } //solve for chemistry
        } //for subCycle

    }
}

void Foam::porousReactionsModels::phreeqcRMModels::updateFluidCompositionOnlyReactions()
{
    Info << "Using SNIA (seq. non-iterative) for transport and chemistry coupling." << nl;

    const volTensorField &Deff = effectiveDispersionTensor();
    word divPhiYiScheme = "div(phi,Yi)";

    std::vector<double> c;
    c.resize(Y_.size()*nxyz_);

    forAll(Y_,i)
    {
            forAll(Y_[i],cellI)
            {
                c[i*nxyz_+cellI]=Y_[i][cellI];
            } //copy flow information to phreeqc concentration
    } //solve transport for the species    i
    updateKineticsParameters();
    Info << "run Phreeqc ... " << nl;
    status = phreeqc_.SetConcentrations(c);
    status = phreeqc_.SetTimeStep(mesh_.time().deltaT().value());
    status = phreeqc_.SetTime(mesh_.time().value());
    Info << "I am going to runcell" << nl;
    status = phreeqc_.RunCells();
    Info << "I finished runcell" << nl;
    status = phreeqc_.GetConcentrations(c);

    forAll(Y_,s)
    {
        forAll(Y_[s],cellI)
        {
            Y_[s][cellI]=c[nxyz_*s+cellI];
        }
    }
    phreeqc_.CloseFiles();
    //if (strangCounter_==3){
    //updateMineralDistribution();
    //updatePorosity();
    updateKineticsParameters();//}
    updatepH();
}
// -------------------------------------------------------------------------//
void Foam::porousReactionsModels::phreeqcRMModels::updateFluidCompositionFromFlow(PtrList<volScalarField>& Y_primary)
{
    Info << "Using SNIA (seq. non-iterative) for transport and chemistry coupling." << nl;

    // const volTensorField &Deff = effectiveDispersionTensor();
    // word divPhiYiScheme = "div(phi,Yi)";

    std::vector<double> c;
    c.resize(Y_primary.size()*nxyz_);

    forAll(Y_primary,i)
    {
            forAll(Y_primary[i],cellI)
            {
                c[i*nxyz_+cellI]=Y_primary[i][cellI];
            } //copy flow information to phreeqc concentration
    } //solve transport for the species    i
    updateKineticsParameters();
    Info << "run Phreeqc ... " << nl;
    status = phreeqc_.SetConcentrations(c);
    status = phreeqc_.SetTimeStep(mesh_.time().deltaT().value());
    status = phreeqc_.SetTime(mesh_.time().value());
    Info << "I am going to runcell" << nl;
    status = phreeqc_.RunCells();
    Info << "I finished runcell" << nl;
    status = phreeqc_.GetConcentrations(c);

    forAll(Y_primary,s)
    {
        forAll(Y_primary[s],cellI)
        {
            Y_primary[s][cellI]=c[nxyz_*s+cellI];
        }
    }
    phreeqc_.CloseFiles();
    //if (strangCounter_==3){
    //updateMineralDistribution();
    //updatePorosity();
    updateKineticsParameters();//}
    updatepH();
}


void Foam::porousReactionsModels::phreeqcRMModels::updateMineralDistribution()
{
//    Info << "Update minerals distribution .... ";
    {
        int n_user = 1;  // equilibrium phases

        status = phreeqc_.SetCurrentSelectedOutputUserNumber(n_user);
        std::vector<double> so;
        status = phreeqc_.GetSelectedOutput(so);

        for (int i = 0; i < phreeqc_.GetSelectedOutputRowCount(); i++)
        {
            forAll(mineralList_,s)
            {
                Ys_[s].storeOldTime();
                if(activatePhaseEquilibrium_[s] == true)
                {
                    Ys_[s][i] = so[2*s*nxyz_ + i]*Vm_[s].value()*1e3;
                }
            }
        }

        n_user = 2;  // kinetics

        status = phreeqc_.SetCurrentSelectedOutputUserNumber(n_user);
      //  std::vector<double> so;
        status = phreeqc_.GetSelectedOutput(so);

        for (int i = 0; i < phreeqc_.GetSelectedOutputRowCount(); i++)
        {
            forAll(mineralList_,s)
            {
                if(activatePhaseEquilibrium_[s] == false)
                {
                    Ys_[s][i] = so[2*s*nxyz_ + i]*Vm_[s].value()*1e3;
                }
            }
        }
    }
    updatePorosityPhreeqc();
//    Info<<"Ok" << endl;
}

// -------------------------------------------------------------------------//

void Foam::porousReactionsModels::phreeqcRMModels::updatepH()
{
//    Info << "Update pH .... ";

    int n_user = 0;

    status = phreeqc_.SetCurrentSelectedOutputUserNumber(n_user);
    std::vector<double> so;
    status = phreeqc_.GetSelectedOutput(so);

    for (int i = 0; i < phreeqc_.GetSelectedOutputRowCount(); i++)
    {
        pH_[i] = so[0.*nxyz_ + i];
    }

    pH_.correctBoundaryConditions();

//    Info<<"Ok" << endl;
}

// -------------------------------------------------------------------------//

void Foam::porousReactionsModels::phreeqcRMModels::activateSelectedOutput()
{
    Info << "Activate SELECTED_OUTPUT ..." ;
    status = phreeqc_.SetSelectedOutputOn(true);

    std::string input;

    // to ouput the pH value
    input =
    " SELECTED_OUTPUT 0; \
        -reset 		false; \
        -pH       true;  \
      END";
    status = phreeqc_.RunString(true, true, true, input.c_str());

    // to ouput the minerals distribution with equilibrium phases
    input =
    " SELECTED_OUTPUT 1; \
        -reset 		false; \
        -equi \
    ";

    forAll(mineralList_,s)
    {
      input +=  mineralList_[s] + " " ;
    }

    input += "; END";

//    Info << "input = " << input << nl << endl;

    status = phreeqc_.RunString(true, true, true, input.c_str());

    // to ouput the minerals distribution with kinetics phases
    input =
    " SELECTED_OUTPUT 2; \
        -reset 		false; \
        -kin \
    ";

    forAll(mineralList_,s)
    {
        input +=  mineralList_[s] + " " ;
    }

    input += "; END";

//    Info << "input = " << input << nl << endl;

    status = phreeqc_.RunString(true, true, true, input.c_str());

    Info << "OK" << nl << endl;
}

// -------------------------------------------------------------------------//
