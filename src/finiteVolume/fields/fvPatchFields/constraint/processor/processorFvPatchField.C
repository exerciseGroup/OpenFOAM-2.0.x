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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p))
{
	//add by Xiaowei:begin
	const myLabelList& faceCells = this->patch().faceCells();
	Time::setFaceCells(&faceCells);
	//add by Xiaowei:end

}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p))
{

    //add by Xiaowei:begin
    const myLabelList& faceCells = this->patch().faceCells();
    Time::setFaceCells(&faceCells);
    //add by Xiaowei:end

}


// Construct by mapping given processorFvPatchField<Type>
template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

	//add by Xiaowei:begin

	const myLabelList& faceCells = this->patch().faceCells();
	Time::setFaceCells(&faceCells);
    //add by Xiaowei:end
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }


	    //add by Xiaowei:begin
	const myLabelList& faceCells = this->patch().faceCells();
	Time::setFaceCells(&faceCells);
    //add by Xiaowei:end
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{
		//add by Xiaowei:begin

		const myLabelList& faceCells = this->patch().faceCells();
		Time::setFaceCells(&faceCells);
		//add by Xiaowei:end
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{
		//add by Xiaowei:begin
		
		const myLabelList& faceCells = this->patch().faceCells();
		Time::setFaceCells(&faceCells);
		//add by Xiaowei:end

}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::~processorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > processorFvPatchField<Type>::patchNeighbourField() const
{
    return *this;
}


template<class Type>
void processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {

	//add by Xiaow: begin
    Foam::Time::enterSec("initEvaluate");
    //add by Xiaow: end
    
    procPatch_.compressedSend(commsType, this->patchInternalField()());

   	//add by Xiaow: begin
    Foam::Time::leaveSec();
    //add by Xiaow: end

	}

}


template<class Type>
void processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        procPatch_.compressedReceive<Type>(commsType, *this);

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }

    }

	
}


template<class Type>
tmp<Field<Type> > processorFvPatchField<Type>::snGrad() const
{
    
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{

    //add by Xiaow: begin
    Foam::Time::enterSec("initInterfaceMatrixUpdate");
    //add by Xiaow: end

    procPatch_.compressedSend
    (
        commsType,
        this->patch().patchInternalField(psiInternal)()
    );

	//add by Xiaow: begin
    Foam::Time::leaveSec();
    //add by Xiaow: end



}


template<class Type>
void processorFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    scalarField pnf
    (
        procPatch_.compressedReceive<scalar>(commsType, this->size())()
    );

    // Transform according to the transformation tensor
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result

    const myLabelList& faceCells = this->patch().faceCells();

//add by Xiaowei:begin
//	Time::setFaceCells(&faceCells);

//add by Xiaowei:end

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
