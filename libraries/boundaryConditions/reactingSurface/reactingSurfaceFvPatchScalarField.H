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
    Foam::reactingSurfaceFvPatchScalarField

Description
    reacting surface boundary condition

SourceFiles
    reactingSurfaceFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef reactingSurfaceFvPatchScalarField_H
#define reactingSurfaceFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class reactingSurfaceFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class reactingSurfaceFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

	word kineticPhase_;

public:

    //- Runtime type information
    TypeName("reactingSurface");


    // Constructors

        //- Construct from patch and internal field
        reactingSurfaceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        reactingSurfaceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given reactingSurfaceFvPatchScalarField
        //  onto a new patch
        reactingSurfaceFvPatchScalarField
        (
            const reactingSurfaceFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        reactingSurfaceFvPatchScalarField
        (
            const reactingSurfaceFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new reactingSurfaceFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        reactingSurfaceFvPatchScalarField
        (
            const reactingSurfaceFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new reactingSurfaceFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

		//get surface master list
		const word& get_kinetic_phase() const
		{
			return kineticPhase_;
		}

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
