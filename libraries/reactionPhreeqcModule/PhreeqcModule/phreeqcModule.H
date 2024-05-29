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

Class
    Foam::phreeqcModule

Description
    class for phreeqc reaction module

SourceFiles
    phreeqcModule.cpp

\*---------------------------------------------------------------------------*/

#ifndef phreeqcModule_H
#define phreeqcModule_H

#include "reactionModule.H"
#include "PhreeqcRM.h"

namespace Foam
{



class phreeqcModule : public reactionModule
{
   

public: 
	phreeqcModule
	(
		const speciesTable&,
		const speciesTable&,
		const speciesTable&,
		const speciesTable&,
		const wordList&,
		const fvMesh& mesh,
		PtrList<volScalarField>&,
		PtrList<volScalarField>&,
		PtrList<volScalarField>&,
		PtrList<volScalarField>&,
		volScalarField&,
		volScalarField&,
		volScalarField&,
        volScalarField&
	);
	virtual ~phreeqcModule();

	//solve reaction step function
	void reactionStep(dimensionedScalar deltaT);

protected:
	//initialise phreeqc module 
	void initialise();

	//index of phreeqc module
	int id_;

	//water saturation
	double* saturation_;

	//index mapping from solution species to phreeqc component
	int* componentSolutionIndex_;
	//index mapping from surface species to phreeqc component
	int* componentSurfaceIndex_;
	//concentration
	double* concentration_;
	//surface concentration
	double* surfConcentration_;
	//surface area
	double* surfArea_;
	//surface potential
	double* surfPotential_;
	//kinetic phase moles
	double* kineticPhaseMoles_;
	//kinetic moles
	double* kineticMoles_;

	//selected output
	double* selectedOutput_;

	//index mappring from selected output to phreeqc selected output
	int* selectedOutputIndex_;

    PhreeqcRM phreeqc_;

    IRM_RESULT status;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
