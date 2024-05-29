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
    Foam::reaction Module

Description
    virtual class for reaction module

SourceFiles
    reactionModule.cpp

\*---------------------------------------------------------------------------*/

#ifndef reactionModule_H
#define reactionModule_H

#include "volFields.H" 
#include "surfaceFields.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "IOdictionary.H"

namespace Foam
{



class reactionModule
{
   

public: 
	reactionModule
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
	virtual ~reactionModule();


	//virtual solve reaction step function
	virtual void reactionStep(dimensionedScalar deltaT);

protected:

	//List of solution species
	speciesTable solutionSpecies_;
	//List of surface species
	speciesTable surfaceSpecies_;
	//List of surface master
	speciesTable surfaceMasters_;
	//List of kinetic phases
	speciesTable kineticPhases_;
	//List of selected output
	wordList selectedOutputNames_;
	//OpenFOAM mesh
	const fvMesh& mesh_;
	//solution species concentration
	PtrList<volScalarField>& Y_;
	//solution species surface concentration
	PtrList<volScalarField>& sY_;
	//kinetic reaction rate
	PtrList<volScalarField>& R_;
	//selected output
	PtrList<volScalarField>& sOut_;
	//water volume fraction
	volScalarField& alpha_;
	//ionic strength
	volScalarField& I_;
	//surface 
	volScalarField& Surf_;
	//surface potential
	volScalarField& psi_;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
