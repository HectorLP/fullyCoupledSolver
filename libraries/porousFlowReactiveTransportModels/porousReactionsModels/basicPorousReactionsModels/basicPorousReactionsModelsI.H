/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
Class
    Foam::basicPorousReactionsModels
Description
    basicPorousReactionsModels module
SourceFiles
    basicPorousReactionsModels.C

    This is the inline function file.
Authors
    17/02/2020 : C. Soulaine
\*---------------------------------------------------------------------------*/

inline Foam::PtrList<Foam::volScalarField>& 
Foam::basicPorousReactionsModels::Y()
{
    return Y_;
}

inline const Foam::PtrList<Foam::volScalarField>& 
Foam::basicPorousReactionsModels::Y() const 
{
    return Y_;
}

inline Foam::volScalarField& Foam::basicPorousReactionsModels::Y(const label i)
{
    return Y_[i];
}

inline const Foam::volScalarField& Foam::basicPorousReactionsModels::Y
(
    const label i
) const
{
    return Y_[i];
}

inline Foam::volScalarField& Foam::basicPorousReactionsModels::Y
(
    const word& specieName
)
{
    return Y_[componentList_[specieName]];
}

inline const Foam::volScalarField& Foam::basicPorousReactionsModels::Y
(
    const word& specieName
) const
{
    return Y_[componentList_[specieName]];
}