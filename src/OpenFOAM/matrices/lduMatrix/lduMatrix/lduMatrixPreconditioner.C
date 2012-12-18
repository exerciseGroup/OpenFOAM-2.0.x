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
#include "diagonalSolver.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //defineRunTimeSelectionTable(lduMatrix::preconditioner, symMatrix);
    //defineRunTimeSelectionTable(lduMatrix::preconditioner, asymMatrix);
    
    //add by NUDT Exercise Group-Xiaowei:begin
    defineRunTimeSelectionTable2(lduMatrix::preconditioner, symMatrix);
    defineRunTimeSelectionTable2(lduMatrix::preconditioner, asymMatrix);
    //add by NUDT Exercise Group-Xiaowei:end
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::lduMatrix::preconditioner::getName
(
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::lduMatrix::preconditioner>
Foam::lduMatrix::preconditioner::New
(
    const solver& sol,
    const dictionary& solverControls
)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = solverControls.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (sol.matrix().symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(name);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::preconditioner::New"
                "(const solver&, const dictionary&)",
                controls
            )   << "Unknown symmetric matrix preconditioner "
                << name << nl << nl
                << "Valid symmetric matrix preconditioners :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            constructorIter()
            (
                sol,
                controls
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(name);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::preconditioner::New"
                "(const solver&, const dictionary&)",
                controls
            )   << "Unknown asymmetric matrix preconditioner "
                << name << nl << nl
                << "Valid asymmetric matrix preconditioners :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            constructorIter()
            (
                sol,
                controls
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::preconditioner::New"
            "(const solver&, const dictionary&)",
            controls
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::preconditioner>(NULL);
    }
}

//add by NUDT Exercise Group-Xiaowei:begin

Foam::autoPtr<Foam::lduPreconditioner>
Foam::lduPreconditioner::New
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
{
    word preconName;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> preconName;
    }
    else
    {
        e.stream() >> preconName;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.symmetric())
    {
        symMatrixConstructorTable1::iterator constructorIter =
            symMatrixConstructorTablePtr_1->find(preconName);

        if (constructorIter == symMatrixConstructorTablePtr_1->end())
        { 
           FatalIOErrorIn
            (
                "lduPreconditioner::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown symmetric matrix preconditioner "
                << preconName << endl << endl
                << "Valid symmetric matrix preconditioners are :" << endl
                << symMatrixConstructorTablePtr_1->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduPreconditioner>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                controls
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable1::iterator constructorIter =
            asymMatrixConstructorTablePtr_1->find(preconName);

        if (constructorIter == asymMatrixConstructorTablePtr_1->end())
        {
            FatalIOErrorIn
            (
                "lduPreconditioner::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown asymmetric matrix preconditioner "
                << preconName << endl << endl
                << "Valid asymmetric matrix preconditioners are :" << endl
                << asymMatrixConstructorTablePtr_1->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduPreconditioner>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                controls
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduPreconditioner::New\n"
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

        return autoPtr<lduPreconditioner>(NULL);
    }
}


           Foam::lduPreconditioner:: preconditioner
            (
                const lduMatrix& matrix,
                const FieldField<Field, scalar>& coupleBouCoeffs,
                const FieldField<Field, scalar>& coupleIntCoeffs,
                const lduInterfaceFieldPtrsList& interfaces
            )
            :
           solver_(*(new diagonalSolver
                        (
                         "abcd",
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                dictionary::null
                        )
                )),
           matrix_(matrix),
           coupleBouCoeffs_(coupleBouCoeffs),
           coupleIntCoeffs_(coupleIntCoeffs),
           interfaces_(interfaces)
           {}
//add by NUDT Exercise Group-Xiaowei:end
// ************************************************************************* //
