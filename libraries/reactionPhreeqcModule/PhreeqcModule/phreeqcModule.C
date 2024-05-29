/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the Carbonate Reservoir Group at Heriot-Watt University. Please visit our
    website for more information <https://carbonates.hw.ac.uk>.
\*---------------------------------------------------------------------------*/


#include "phreeqcModule.H"
#include "RM_interface_C.h"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"
#include "reactingWallFvPatchScalarField.H"
#include "reactingSurfaceFvPatchScalarField.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::phreeqcModule::initialise()
{
	//get number of cells
    Info << "Start to initialize the reaction module." << endl;
	int ncells = mesh_.cells().size();
    Info << "Cell number is " << ncells << endl;
	// id_ = RM_Create(ncells, 1);
    phreeqc_(ncells, 1);
    // int tmpId = id_;
    bool workers = true;
    bool initial_phreeqc = true;
    bool utility = true;
    Info << "The number for seting is " << tmpId << endl;
	//set porosity=1.0
	// double* porosity = new double[ncells];
	// for (int i = 0; i < ncells; i++) porosity[i] = 1.0;
	// RM_SetPorosity(id_, porosity);
    std::vector<double> por;
    por.resize(ncells)
    status = phreeqc_.SetPorosity(por);
	// delete[] porosity;

	Info << "Set the solution concentration unit." << endl;
	//concentration=mol/L
	// RM_SetUnitsSolution(id_, 2);
    status = phreeqc_.SetUnitsSolution(2);

	//Save species for multi-species transport
    Info << "Set the species for saving." << endl;
	// RM_SetSpeciesSaveOn(id_, true);
    status = phreeqc_.SetSpeciesSaveOn(true);

	//Save selected output
    Info << "Set the select output part." << endl;
	// RM_SetSelectedOutputOn(id_,true);
    status = phreeqc_.SetSelectedOutputOn(true);

	//load phreeqc database
    Info << "Load the database." << endl;
	// RM_LoadDatabase(id_, "constant/GeoChem.dat");
    status = phreeqc_.LoadDatabase("constant/GeoChem.dat");

    Info << "Load the Input file." << endl;
    status = phreeqc_.RunFile(workers, initial_phreeqc, utility, "constant/phreeqcReactions");

    // RM_RunFile(id_, 1, 1, 0, "constant/phreeqcReactions");

	// RM_LoadDatabase(id_, "constant/GeoChem.dat");

    // RM_RunFile(id_, 1, 1, 1, "constant/phreeqcReactions");
	//load reaction
	

	//load reaction


	//solution component list
	// std::ostringstream oss;
    std::string input;
    // Info <<  oss.str().c_str() << endl;

	//init surface and solution
	forAll(mesh_.cells(),celli)
	{
		// oss << "COPY solution 0 " << celli  << "\n";
		// oss << "END" << "\n";
        input += "COPY solution 0 " +  std::to_string(celli) + "\n";
        input += "END \n";
    }

	//equilibrate surface with solution
	forAll(mesh_.cells(),celli)
	{
		// oss << "SURFACE " << celli << "\n";
		// oss << "-equilibrate with solution " << celli << "\n";
        input += "SURFACE " + std::to_string(celli) + "\n";
        input += "-equilibrate with solution " + std::to_string(celli) + "\n";
		forAll(surfaceMasters_, j)
		{
			double area = 0.0;
			double mole = 0.0;
	        const volScalarField::Boundary& Surfbf = Surf_.boundaryField();
			forAll(Surfbf,patchi)
			{
				if (Surfbf[patchi].type() == "reactingWall")
				{
					const reactingWallFvPatchScalarField& Surfcap = refCast<const reactingWallFvPatchScalarField>(Surfbf[patchi]);
					const labelList& cellOwner = Surfcap.patch().faceCells();
					const surfaceScalarField& magSf = mesh_.magSf();
					const wordList& masters = Surfcap.get_surface_masters();
					const scalarList& density   = Surfcap.get_density();//mol/m^2
					forAll(masters,i)
					{
						if (masters[i]==surfaceMasters_[j])
						{
							forAll(Surfbf[patchi],facei)
							{
								if (cellOwner[facei]==celli)
								{
									mole+=density[i]*magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//mol/L
									area+=magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//m^2/L
								}
							}
						}
					}
				}
			}	
            // if (area>0) oss << surfaceMasters_[j] << "  " << mole << " " << area  << " 1" << "\n";
            // else oss << surfaceMasters_[j] << "  " << "0 1 1" << "\n";
            if (area>0) input += surfaceMasters_[j] + "  " + std::to_string(mole) + " " + std::to_string(area) + " " + std::to_string(1) + "\n";
            else input += surfaceMasters_[j] + "  " + "0 1 1" + "\n"; 
		}
		// oss << "END" << "\n";
        input += "END \n"; 
	}

	//equilibrate surface with solution
	// forAll(mesh_.cells(),celli)
	// {
	// 	oss << "KINETICS " << celli << "\n";
	// 	forAll(kineticPhases_,j)
	// 	{
	// 		double area = 0.0;
	// 		double mole = 0.0;
	//         const volScalarField::Boundary& Surfbf = Surf_.boundaryField();
	// 		forAll(Surfbf,patchi)
	// 		{
	// 			if (Surfbf[patchi].type() == "reactingSurface")
	// 			{
	// 				const reactingSurfaceFvPatchScalarField& Surfcap = refCast<const reactingSurfaceFvPatchScalarField>(Surfbf[patchi]);
	// 				const labelList& cellOwner = Surfcap.patch().faceCells();
	// 				const surfaceScalarField& magSf = mesh_.magSf();
	// 				word phase = Surfcap.get_kinetic_phase();
	// 				if (phase==kineticPhases_[j])
	// 				{
	// 					forAll(Surfbf[patchi],facei)
	// 					{
	// 						if (cellOwner[facei]==celli)
	// 						{
	// 							mole+=1.0*magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//mol/L
	// 							area+=magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//m^2/L
	// 						}
	// 					}
	// 				}			
	// 			}
	// 		}				
		
	// 		if (mole>0)
	// 		{ 
	// 			oss << kineticPhases_[j] << "\n";
	// 		    oss << "-tol 1e-8" << "\n";
	// 			oss << "-m0 " << mole << "\n";
	// 			oss << "-m " << mole << "\n";
	// 			oss << "-parms " << area/mole*1e4 << "\n";//convert m^2/mol to cm^2/mol
	// 		}	
	// 	}
    //     oss << "END" << "\n";
    //     if (celli == 559999)
    //     {
    //         oss << "This loop should be ended." << "\n";
    //     }
	// }

    Info << "Load the reaction into PhreeqcRM." << endl;
    // Info <<  oss.str().c_str() << endl;
    // RM_RunString(id_, 1, 1, 0, oss.str().c_str());
    status = phreeqc_.RunString(true, true, true, input.c_str());
    
    // Info << "Show phreeqc keywords." <<endl;
	//display phreeqc keywords

	//run phreeqc keywords
	
    Info << "Initialize Phreeqc module." << endl;
	//init Phreeqc worker module
	// int* ic1 = (int *)malloc((size_t)(7 * ncells * sizeof(int)));
    std::vector<int> ic1;
    ic1.resize(ncells, -1);
    Info << "Initialize Phreeqc cell by cell." << endl;
	forAll(mesh_.cells(),celli)
	{
		ic1[celli] = celli;               // Solution  i
		ic1[ncells + celli] = -1;      // Equilibrium phases none
		ic1[2 * ncells + celli] = -1;       // Exchange none
		ic1[3 * ncells + celli] = celli;      // Surface i
		ic1[4 * ncells + celli] = -1;      // Gas phase none
		ic1[5 * ncells + celli] = -1;      // Solid solutions none
		ic1[6 * ncells + celli] = -1;      // Kinetics i
	}
    Info << "Initialize using phreeqc2module." << endl;
	// RM_InitialPhreeqc2Module(id_, ic1, 0, 0);
    status = phreeqc_.InitialPhreeqc2Module(ic1);
	// free(ic1);

	//Run Phreeqc to init concentration
	// RM_RunCells(id_);
    status = phreeqc_.RunCells();

	//find components
	// RM_FindComponents(id_);
    status = phreeqc.FindComponents();
	
	//get number of solution species
	// int nsol = RM_GetSpeciesCount(id_);
    int nspecies = phreeqc_.GetSpeciesCount();

	//display number of silution species
	Info << "number of solution:" << nsol << "\n";

	//set solution speciesvector
	// char** components = (char **)malloc((size_t)(nsol * sizeof(char *)));

	//get solution species name
	// for (int i = 0; i < nsol; i++)
	// {
	// 	components[i] = (char *)malloc((size_t)(20 * sizeof(char *)));
	// 	RM_GetSpeciesName(id_, i, components[i], 20);
	// }
    const std::vector<std::string> &species = phreeqc_.GetSpeciesNames();

	//get number of surface species
	// int nsurf = RM_GetSurfaceSpeciesCount(id_);
    int nsurf = phreeqc_.GetSurfaceSpeciesCount();

	//display number of surface species
	Info << "number of surface:" << nsurf << "\n";

	//set solution speciesvector
	// char** surfComponents = (char **)malloc((size_t)(nsurf * sizeof(char *)));

	// //get solution species name
	// for (int i = 0; i < nsurf; i++)
	// {
	// 	surfComponents[i] = (char *)malloc((size_t)(20 * sizeof(char *)));
	// 	RM_GetSurfaceSpeciesName(id_, i, surfComponents[i], 20);
	// }
    const std::vector<std::string> surfComponents = phreeqc_.GetSurfaceSpcies();
	
	//save component index map
    Info << "Save the component index." << endl;
	componentSolutionIndex_ = (int*)malloc((size_t)(solutionSpecies_.size() * sizeof(int)));
	componentSurfaceIndex_  = (int*)malloc((size_t)(surfaceSpecies_.size() * sizeof(int)));
	// forAll(solutionSpecies_, i)
	// {
	// 	for (int j = 0; j < nsol; j++)
	// 	{
	// 		std::string component = components[j];
	// 		if (component == solutionSpecies_[i])
	// 		{
	// 			componentSolutionIndex_[i] = j;
	// 		}
	// 	}
	// }
    forAll(solutionSpecies_, i)
    {
        for (int j = 0; j < nspecies; ++j)
        {
            std::string component = species[j];
            if (component == solutionSpecies_[i])
            {
                componentSolutionIndex_[i] = j;
            }
        }
    }
    Info << "Save the surface species index." << endl;
	forAll(surfaceSpecies_, i)
	{
		for (int j = 0; j < nsurf; j++)
		{
			std::string component = surfComponents[j];
			if (component == surfaceSpecies_[i])
			{
				componentSurfaceIndex_[i] = j;
			}
		}
	}

	//concentration, adsorption and surface potential
	if (nsol>0) concentration_      = (double *)malloc((size_t)(nsol * ncells * sizeof(double)));
	if (nsurf>0) 
	{
        Info << "Initialize the size of array for surface." << endl;
		surfConcentration_  = (double *)malloc((size_t)(nsurf * ncells * sizeof(double)));
		surfArea_           = (double *)malloc((size_t)(ncells * sizeof(double)));
		surfPotential_      = (double *)malloc((size_t)(ncells * sizeof(double)));
	}
    std::vector<double>& concentrationSpecies;
    
    std:vector<double>& surfConcentration;
    
    std::vector<double>& surfPotential;
    
    std::vector<double>& surfArea;
     
	if (nspecies > 0)
    {
        concentrationSpecies.resize(nspecies * ncells);
    }
    if (nsurf > 0)
    {
        surfConcentration_.resize(nsurf * ncells);
        surfPotential_.resize(ncells);
        surfArea_.resize(ncells);
    }
	// RM_GetSpeciesConcentrations(id_, concentration_);
    status = phreeqc_.GetSpeciesConcentration(concentrationSpecies);
    status = phreeqc_.GetSurfaceSpeciesConcentrations(surfConcentration_);
    status = phreeqc_.GetSurfacePotential(surfPotential);
    status = phreeqc_.GetSurfaceArea(surfaceArea);

	// RM_GetSurfaceSpeciesConcentrations(id_, surfConcentration_);
	// RM_GetSurfaceArea(id_,"Surf",surfArea_);
	// RM_GetSurfacePotential(id_,"Surf",surfPotential_);

	//set water saturation for Phreeqc module
    Info << "Initilize the chemical concentration." << endl;
	// saturation_ = (double *)malloc((size_t)(ncells * sizeof(double)));
	// forAll(mesh_.cells(), celli)
	// {
	// 	if (alpha_[celli] >1e-3) saturation_[celli] = alpha_[celli];
	// 	else  saturation_[celli]=0;
	// }
	// RM_SetSaturation(id_, saturation_);
    std::vector<double>& saturation_;
    saturation_.resize(ncells);
    forAll(mesh_.cells(), celli)
    {
        if (alpha_[celli] >1e-3) saturation_[celli] = alpha_[celli];
        else saturation_[celli] = 0.0;
    }
	//get concentration from OpenFOAM
    Info << "Initialize the chemical concentration from OpenFOAM." << endl;
	forAll(solutionSpecies_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			// concentration_[componentSolutionIndex_[i] * ncells + celli] = Yi[celli]/1000;
            concentrationSpecies[componentSolutionIndex_[i] * ncells + celli] = Yi[celli] / 1000;

		}
	}

	//get surface concentration from OpenFoam
    Info << "Initialize the chemical concentration from OpenFOAM." << endl;
	forAll(surfaceSpecies_, i)
	{
		forAll(mesh_.cells(), celli)
		{
			surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] = 0.0;
		}
		volScalarField& Yi = sY_[i];
	    forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				const scalarField& Yfaces = Yi.boundaryField()[patchi];
				const surfaceScalarField& magSf = mesh_.magSf();
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (surfArea_[cellOwner[facei]]>0)
					{
						surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]] += Yfaces[facei] / surfArea_[cellOwner[facei]] * magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;
					}
				}
			}
		}
	}


	//kinetic phases moles
	// int nphases = kineticPhases_.size();

	// Info << "nphases:" << nphases << endl;
	// if (nphases>0)
	// {
	// 	kineticPhaseMoles_ = (double *)malloc((size_t)(nphases*ncells*sizeof(double)));
	// 	kineticMoles_ = (double *)malloc((size_t)(nphases*ncells*sizeof(double)));
	// }

	// forAll(kineticPhases_,i)
	// {
	// 	RM_GetKineticPhaseMoles(id_,kineticPhases_[i],&kineticPhaseMoles_[i*ncells]);
	// }
 	
	//save selected output
	// int ncol = RM_GetSelectedOutputColumnCount(id_);
    int ncol = phreeqc_.GetSelectedOutputColumnCount();
    selectedOutput_      = (double *)malloc((size_t)(ncol*ncells*sizeof(double)));
	selectedOutputIndex_ = (int*)malloc((size_t)(selectedOutputNames_.size() * sizeof(int)));
	// char* heading = (char*)malloc((size_t)(30*sizeof(char)));
    std::string& heading

	forAll(selectedOutputNames_, i)
	{
		for (int j = 0; j < ncol; j++)
		{
			// RM_GetSelectedOutputHeading(id_,j,heading,30);
            status = phreeqc_.GetSelectedOutputHeading(j, heading);
			std::string head=heading;
			if (head == selectedOutputNames_[i])
			{
				selectedOutputIndex_[i] = j;
			}
			
		}
	}
	// free(heading);

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phreeqcModule::phreeqcModule
(
	const speciesTable& solutionSpecies,
	const speciesTable& surfaceSpecies,
	const speciesTable& surfaceMasters,
	const speciesTable& kineticPhases,
	const wordList& selectedOutputNames,
	const fvMesh& mesh,
	PtrList<volScalarField>& Y,
	PtrList<volScalarField>& sY,
	PtrList<volScalarField>& R,
	PtrList<volScalarField>& sOut,
	volScalarField& alpha,
	volScalarField& I,
	volScalarField& Surf,
	volScalarField& psi
)
:
reactionModule(solutionSpecies, surfaceSpecies, surfaceMasters,kineticPhases,selectedOutputNames,mesh,Y,sY,R,sOut,alpha,I,Surf,psi),
saturation_(NULL),
componentSolutionIndex_(NULL),
componentSurfaceIndex_(NULL),
concentration_(NULL),
surfConcentration_(NULL),
surfArea_(NULL),
surfPotential_(NULL),
kineticPhaseMoles_(NULL),
kineticMoles_(NULL),
selectedOutput_(NULL),
selectedOutputIndex_(NULL)
{
	initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::phreeqcModule::~phreeqcModule
(
)
{
	free(saturation_);
	free(componentSolutionIndex_);
	free(componentSurfaceIndex_);
	free(concentration_);
	free(surfConcentration_);
	free(surfArea_);
	free(surfPotential_);
	free(kineticPhaseMoles_);
	free(kineticMoles_);
	free(selectedOutput_);
	free(selectedOutputIndex_);
}


// * * * * * * * * * * * * * * * * solve reaction step function  * * * * * * * * * * * * * * //
void Foam::phreeqcModule::reactionStep
(
	dimensionedScalar deltaT
)
{
	//get number of cells
	int ncells = mesh_.cells().size();

	//get water saturation after transport
	forAll(mesh_.cells(), celli)
	{
		if (alpha_[celli] >1e-3) saturation_[celli] = alpha_[celli];
		else  saturation_[celli]=0;
	}
	RM_SetSaturation(id_, saturation_);

	//get concentration after transport
	forAll(solutionSpecies_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			concentration_[componentSolutionIndex_[i] * ncells + celli] = Yi[celli]/1000;//mol/L
		}
	}


	//get surface concentration after transport
	forAll(surfaceSpecies_, i)
	{
		forAll(mesh_.cells(), celli)
		{
			surfConcentration_[componentSurfaceIndex_[i] * ncells + celli] = 0.0;
		}
		volScalarField& Yi = sY_[i];
	    forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				const scalarField& Yfaces = Yi.boundaryField()[patchi];
				const surfaceScalarField& magSf = mesh_.magSf();
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (surfArea_[cellOwner[facei]]>0)
					{
						surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]] += Yfaces[facei] / surfArea_[cellOwner[facei]] * magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;
					}
				}
			}
		}
	}

	//equilibrate surface with solution
	// forAll(mesh_.cells(),celli)
	// {
	// 	forAll(kineticPhases_,j)
	// 	{
	// 		double area = 0.0;
	// 		double mole = 0.0;
	//         const volScalarField::Boundary& Surfbf = Surf_.boundaryField();
	// 		forAll(Surfbf,patchi)
	// 		{
	// 			if (Surfbf[patchi].type() == "reactingSurface")
	// 			{
	// 				const reactingSurfaceFvPatchScalarField& Surfcap = refCast<const reactingSurfaceFvPatchScalarField>(Surfbf[patchi]);
	// 				const labelList& cellOwner = Surfcap.patch().faceCells();
	// 				const surfaceScalarField& magSf = mesh_.magSf();
	// 				word phase = Surfcap.get_kinetic_phase();
	// 				if (phase==kineticPhases_[j])
	// 				{
	// 					forAll(Surfbf[patchi],facei)
	// 					{
	// 						if (cellOwner[facei]==celli)
	// 						{
	// 							mole+=1.0*magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000 ;//mol/L
	// 							area+=magSf.boundaryField()[patchi][facei] / mesh_.V()[cellOwner[facei]] / 1000;//m^2/L
	// 						}
	// 					}
	// 				}			
	// 			}
	// 		}				
		
	// 		kineticPhaseMoles_[j*ncells+celli]=mole;
	// 	}
	// }

	//set concentration for Phreeqc module
	RM_SpeciesConcentrations2Module(id_, concentration_);
	RM_SurfaceSpeciesConcentrations2Module(id_, surfConcentration_);
	// forAll(kineticPhases_,i)
	// {
	// 	RM_SetKineticPhaseMoles(id_,kineticPhases_[i],&kineticPhaseMoles_[i*ncells]);
	// }

	RM_SetTimeStep(id_,deltaT.value());

	//Run phreeqc
	RM_RunCells(id_);

	//get concentration after reactions
	RM_GetSpeciesConcentrations(id_, concentration_);
	RM_GetSurfaceSpeciesConcentrations(id_, surfConcentration_);

	//get surface potential
	RM_GetSurfacePotential(id_,"Surf",surfPotential_);

	//get ionic strength
	RM_GetSolutionIonicStrength(id_, &I_[0]);

	//get reaction rate and new vector composition
	forAll(solutionSpecies_, i)
	{
		volScalarField& Yi = Y_[i];
		forAll(mesh_.cells(), celli)
		{
			if (saturation_[celli]>1e-3)
			{
				Yi[celli]  = concentration_[componentSolutionIndex_[i] * ncells + celli]*1000;//mol/m3
			}
		}
	}

	//get surface concentration
	forAll(surfaceSpecies_, i)
	{
		volScalarField& Yi = sY_[i];
		forAll(Yi.boundaryField(), patchi)
		{
			if (Surf_.boundaryField()[patchi].type()=="reactingWall")
			{
				const labelList& cellOwner = Yi.boundaryField()[patchi].patch().faceCells();
				scalarField& Yfaces   = Yi.boundaryFieldRef()[patchi];
				scalarField& psifaces = psi_.boundaryFieldRef()[patchi];
				forAll(Yi.boundaryField()[patchi], facei)
				{
					if (saturation_[cellOwner[facei]]>1e-3)
					{
						Yfaces[facei] = surfConcentration_[componentSurfaceIndex_[i] * ncells + cellOwner[facei]];
						psifaces[facei]   = surfPotential_[cellOwner[facei]];
					}
				}
			}
		}
	}


	// forAll(kineticPhases_,i)
	// {
	// 	volScalarField& Ri = R_[i];
	// 	RM_GetKineticMoles(id_,kineticPhases_[i],&kineticMoles_[i*ncells]);
	// 	forAll(Surf_.boundaryField(), patchi)
	// 	{
	// 		if (Surf_.boundaryField()[patchi].type()=="reactingSurface")
	// 		{
	// 			const labelList& cellOwner = Ri.boundaryField()[patchi].patch().faceCells();
	// 			scalarField& Rfaces   = Ri.boundaryFieldRef()[patchi];
	// 			const surfaceScalarField& magSf = mesh_.magSf();
	// 			forAll(Ri.boundaryField()[patchi], facei)
	// 			{
	// 				if (saturation_[cellOwner[facei]]>1e-3)
	// 				{
	// 					Rfaces[facei] = kineticMoles_[i*ncells+cellOwner[facei]]/deltaT.value()
	// 								*mesh_.V()[cellOwner[facei]]/magSf.boundaryField()[patchi][facei]*1000;//mol/m2/s
	// 				}
	// 			}
	// 		}
	// 	}
		
	// }

	
	//get selected output
	RM_GetSelectedOutput(id_,selectedOutput_);

	//get new selected output vector
	forAll(selectedOutputNames_, i)
	{
		volScalarField& sOuti = sOut_[i];
		forAll(mesh_.cells(), celli)
		{
			if (saturation_[celli]>1e-3)
			{
				sOuti[celli] = selectedOutput_[selectedOutputIndex_[i]*ncells+celli];
			}
		}	
	}


}

// ************************************************************************* //
