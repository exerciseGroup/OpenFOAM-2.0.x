/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
/*
// delete by NUDT Exercise Group-RXG: begin
    defineRunTimeSelectionTable(lduMatrix::smoother, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::smoother, asymMatrix);
// delete by NUDT Exercise Group-RXG: end
*/

// add by NUDT Exercise Group-RXG: begin
    defineRunTimeSelectionTable2(lduMatrix::smoother, symMatrix);
    defineRunTimeSelectionTable2(lduMatrix::smoother, asymMatrix);
// add by NUDT Exercise Group-RXG: end
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word
Foam::lduMatrix::smoother::getName
(
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::lduMatrix::smoother> Foam::lduMatrix::smoother::New
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    // not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(name);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New", solverControls
            )   << "Unknown symmetric matrix smoother "
                << name << nl << nl
                << "Valid symmetric matrix smoothers are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(name);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New", solverControls
            )   << "Unknown asymmetric matrix smoother "
                << name << nl << nl
                << "Valid asymmetric matrix smoothers are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::smoother::New", solverControls
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::smoother>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::smoother::smoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    fieldName_(fieldName),
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaceIntCoeffs_(interfaceIntCoeffs),
    coupleBouCoeffs_(interfaceIntCoeffs),// add by NUDT Exercise Group-RXG
    coupleIntCoeffs_(interfaceIntCoeffs),// add by NUDT Exercise Group-RXG
    interfaces_(interfaces)
{}


// add by NUDT Exercise Group-RXG: begin
Foam::autoPtr<Foam::lduMatrix::smoother> Foam::lduMatrix::smoother::New
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
{
    word smootherName;

    // Handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> smootherName;
    }
    else
    {
        e.stream() >> smootherName;
    }

    // Not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.symmetric())
    {
        symMatrixConstructorTable1::iterator constructorIter =
            symMatrixConstructorTablePtr_1->find(smootherName);

        if (constructorIter == symMatrixConstructorTablePtr_1->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown symmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid symmetric matrix smoothers are :" << endl
                << symMatrixConstructorTablePtr_1->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduSmoother>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable1::iterator constructorIter =
            asymMatrixConstructorTablePtr_1->find(smootherName);

        if (constructorIter == asymMatrixConstructorTablePtr_1->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown asymmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid asymmetric matrix smoothers are :" << endl
                << asymMatrixConstructorTablePtr_1->toc()
    << exit(FatalIOError);
        }

        return autoPtr<lduSmoother>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::smoother::New\n"
            "(\n"
            "    const lduMatrix& matrix,\n"
            "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
            "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
            "    const lduInterfaceFieldPtrsList& interfaces,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduSmoother>(NULL);
    }
}

// add by NUDT Exercise Group-RXG: end

// add by NUDT Exercise Group-RXG: begin
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::lduMatrix::smoother::smoother
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    matrix_(matrix),
    coupleBouCoeffs_(coupleBouCoeffs),
    coupleIntCoeffs_(coupleIntCoeffs),
    interfaceBouCoeffs_(coupleBouCoeffs),// add by NUDT Exercise Group-RXG
    interfaceIntCoeffs_(coupleIntCoeffs),// add by NUDT Exercise Group-RXG
    interfaces_(interfaces)
{}


// ************************************************************************* //

// add by NUDT Exercise Group-RXG: end
// ************************************************************************* //
