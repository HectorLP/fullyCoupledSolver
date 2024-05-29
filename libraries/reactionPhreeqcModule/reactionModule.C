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


#include "reactionModule.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionModule::reactionModule
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
solutionSpecies_(solutionSpecies),
surfaceSpecies_(surfaceSpecies),
surfaceMasters_(surfaceMasters),
kineticPhases_(kineticPhases),
selectedOutputNames_(selectedOutputNames),
mesh_(mesh),
Y_(Y),
sY_(sY),
R_(R),
sOut_(sOut),
alpha_(alpha),
I_(I),
Surf_(Surf),
psi_(psi)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::reactionModule::~reactionModule
(
)
{}


// * * * * * * * * * * * * * * * * solve reaction step function  * * * * * * * * * * * * * * //
void Foam::reactionModule::reactionStep
(
dimensionedScalar deltaT
)
{
	Info << "WARNING: no water phase in system" << endl;
}
// ************************************************************************* //
