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
    defineRunTimeSelectionTable(lduMatrix::solver, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::solver, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduMatrix::solver> Foam::lduMatrix::solver::New
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
{
    const word name(solverControls.lookup("solver"));

    if (matrix.diagonal())
    {
        return autoPtr<lduMatrix::solver>
        (
            new diagonalSolver
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
            )
        );
    }
    else if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(name);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::solver::New", solverControls
            )   << "Unknown symmetric matrix solver " << name << nl << nl
                << "Valid symmetric matrix solvers are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
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
                "lduMatrix::solver::New", solverControls
            )   << "Unknown asymmetric matrix solver " << name << nl << nl
                << "Valid asymmetric matrix solvers are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::solver::New", solverControls
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::solver>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::solver::solver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    fieldName_(fieldName),
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaceIntCoeffs_(interfaceIntCoeffs),
    coupleBouCoeffs_(interfaceBouCoeffs),//add by NUDT Exercise Group-Xiaowei
    coupleIntCoeffs_(interfaceIntCoeffs),//add by NUDT Exercise Group-Xiaowei
    minIter_(0),
    interfaces_(interfaces),
    controlDict_(solverControls)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lduMatrix::solver::readControls()
{
    maxIter_   = controlDict_.lookupOrDefault<label>("maxIter", 1000);
    tolerance_ = controlDict_.lookupOrDefault<scalar>("tolerance", 1e-6);
    relTol_    = controlDict_.lookupOrDefault<scalar>("relTol", 0);
}


void Foam::lduMatrix::solver::read(const dictionary& solverControls)
{
    controlDict_ = solverControls;
    readControls();
}


Foam::scalar Foam::lduMatrix::solver::normFactor
(
    const scalarField& psi,
    const scalarField& source,
    const scalarField& Apsi,
    scalarField& tmpField
) const
{
    // --- Calculate A dot reference value of psi
    matrix_.sumA(tmpField, interfaceBouCoeffs_, interfaces_);
    tmpField *= gAverage(psi);

    return gSum(mag(Apsi - tmpField) + mag(source - tmpField)) + matrix_.small_;

    // At convergence this simpler method is equivalent to the above
    // return 2*gSumMag(source) + matrix_.small_;
}

//add by NUDT Exercise Group-RXG: begin
Foam::scalar Foam::lduMatrix::solver::normFactor
(
    const scalarField& x,
    const scalarField& b,
    const scalarField& Ax,
    scalarField& tmpField,
    const direction cmpt
) const
{
    // Calculate A dot reference value of x
//     matrix_.sumA(tmpField, coupleBouCoeffs_, interfaces_);
//     tmpField *= gAverage(x);

    // Calculate normalisation factor using full multiplication
    // with mean value.  HJ, 5/Nov/2007
    scalar xRef = gAverage(x);
    matrix_.Amul
    (
        tmpField,
        scalarField(x.size(), xRef),
        coupleBouCoeffs_,
        interfaces_,
        cmpt
    );

    return gSum(mag(Ax - tmpField) + mag(b - tmpField)) + matrix_.small_;

    // At convergence this simpler method is equivalent to the above
    // return 2*gSumMag(b) + matrix_.small_;
}
//add by NUDT Exercise Group-RXG: end

//add by NUDT Exercise Group-Xiaowei:begin
Foam::scalar Foam::lduMatrix::solver::normFactor
(
 // const scalarField& x,// delete by NUDT Exercise Group-RXG
    scalarField& x,// add by NUDT Exercise Group-RXG
    const scalarField& b,
    const direction cmpt
) const
{
    scalarField wA(x.size());
    scalarField tmpField(x.size());

    matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);

    return normFactor(x, b, wA, tmpField, cmpt);
}

//add by NUDT Exercise Group-Xiaowei:end

// add by NUDT Exercise Group-RXG: begin
bool Foam::lduMatrix::solver::stop
(
    lduMatrix::solverPerformance& solverPerf
) const
{
    if (solverPerf.nIterations() < minIter_)
    {
        return false;
    }

    if
    (
        solverPerf.nIterations() >= maxIter_
     || solverPerf.checkConvergence(tolerance_, relTol_)//modified by NUDT Exercise Group-Xiaowei
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

// add by NUDT Exercise Group-RXG: end


// ************************************************************************* //
