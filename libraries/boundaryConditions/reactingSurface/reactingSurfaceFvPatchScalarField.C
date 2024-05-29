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

#include "reactingSurfaceFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactingSurfaceFvPatchScalarField::reactingSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::reactingSurfaceFvPatchScalarField::reactingSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
	kineticPhase_(dict.lookup("kineticPhase"))
{
}


Foam::reactingSurfaceFvPatchScalarField::reactingSurfaceFvPatchScalarField
(
    const reactingSurfaceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
	kineticPhase_(ptf.kineticPhase_)
{
}


Foam::reactingSurfaceFvPatchScalarField::reactingSurfaceFvPatchScalarField
(
    const reactingSurfaceFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
	kineticPhase_(ptf.kineticPhase_)
{}


Foam::reactingSurfaceFvPatchScalarField::reactingSurfaceFvPatchScalarField
(
    const reactingSurfaceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
	kineticPhase_(ptf.kineticPhase_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactingSurfaceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
	os.writeKeyword("kineticPhase") << kineticPhase_ << token::END_STATEMENT<< endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        reactingSurfaceFvPatchScalarField
    );
}

// ************************************************************************* //
