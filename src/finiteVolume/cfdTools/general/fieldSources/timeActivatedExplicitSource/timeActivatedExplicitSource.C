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

#include "TimeActivatedExplicitSource.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::TimeActivatedExplicitSource<Type>::
selectionModeTypeNames_
(
    IStringStream("(points cellSet cellZone all)")()
);


template<class Type>
const Foam::wordList Foam::TimeActivatedExplicitSource<Type>::
volumeModeTypeNames_
(
    IStringStream("(absolute specific)")()
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
typename Foam::TimeActivatedExplicitSource<Type>::selectionModeType
Foam::TimeActivatedExplicitSource<Type>::wordToSelectionModeType
(
    const word& smtName
) const
{
    forAll(selectionModeTypeNames_, i)
    {
        if (smtName == selectionModeTypeNames_[i])
        {
            return selectionModeType(i);
        }
    }

    FatalErrorIn
    (
        "TimeActivatedExplicitSource<Type>::selectionModeType"
        "TimeActivatedExplicitSource<Type>::wordToSelectionModeType"
        "("
            "const word&"
        ")"
    )   << "Unknown selectionMode type " << smtName
        << ". Valid selectionMode types are:" << nl << selectionModeTypeNames_
        << exit(FatalError);

    return selectionModeType(0);
}


template<class Type>
typename Foam::TimeActivatedExplicitSource<Type>::volumeModeType
Foam::TimeActivatedExplicitSource<Type>::wordToVolumeModeType
(
    const word& vmtName
) const
{
    forAll(volumeModeTypeNames_, i)
    {
        if (vmtName == volumeModeTypeNames_[i])
        {
            return volumeModeType(i);
        }
    }

    FatalErrorIn
    (
        "TimeActivatedExplicitSource<Type>::volumeModeType"
        "TimeActivatedExplicitSource<Type>::wordToVolumeModeType(const word&)"
    )   << "Unknown volumeMode type " << vmtName
        << ". Valid volumeMode types are:" << nl << volumeModeTypeNames_
        << exit(FatalError);

    return volumeModeType(0);
}


template<class Type>
Foam::word Foam::TimeActivatedExplicitSource<Type>::selectionModeTypeToWord
(
    const selectionModeType& smtType
) const
{
    if (smtType > selectionModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return selectionModeTypeNames_[smtType];
    }
}


template<class Type>
Foam::word Foam::TimeActivatedExplicitSource<Type>::volumeModeTypeToWord
(
    const volumeModeType& vmtType
) const
{
    if (vmtType > volumeModeTypeNames_.size())
    {
        return "UNKNOWN";
    }
    else
    {
        return volumeModeTypeNames_[vmtType];
    }
}


template<class Type>
void Foam::TimeActivatedExplicitSource<Type>::setSelection
(
    const dictionary& dict
)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "TimeActivatedExplicitSource::setSelection(const dictionary&)"
            )   << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::TimeActivatedExplicitSource<Type>::setFieldData
(
    const dictionary& dict,
    const wordList& fieldNames
)
{
    dict.lookup("fieldData") >> fieldData_;
    labelList localFieldIds(fieldData_.size(), -1);
    forAll(fieldNames, i)
    {
        forAll(fieldData_, j)
        {
            const word& fdName = fieldData_[j].first();
            if (fdName == fieldNames[i])
            {
                fieldIds_[i] = j;
                localFieldIds[j] = i;
                break;
            }
        }
    }
    forAll(localFieldIds, i)
    {
        if (localFieldIds[i] < 0)
        {
            FatalErrorIn
            (
                "TimeActivatedExplicitSource<Type>::setFieldData"
                "("
                    "const dictionary&, "
                    "const wordList&"
                ")"
            )   << "Field " << fieldData_[i].first() << " not found in "
                << "field list. Available fields are: " << nl << fieldNames
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::TimeActivatedExplicitSource<Type>::setCellSet()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl;
    switch (selectionMode_)
    {
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label cellI = mesh_.findCell(points_[i]);
                if (cellI >= 0)
                {
                    selectedCells.insert(cellI);
                }

                label globalCellI = returnReduce(cellI, maxOp<label>());
                if (globalCellI < 0)
                {
                    WarningIn("TimeActivatedExplicitSource<Type>::setCellIds()")
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;
                }
            }

            cells_ = selectedCells.toc();

            break;
        }
        case smCellSet:
        {
            Info<< indent << "- selecting cells using cellSet "
                << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent << "- selecting cells using cellZone "
                << cellSetName_ << endl;
            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorIn("TimeActivatedExplicitSource<Type>::setCellIds()")
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorIn("TimeActivatedExplicitSource<Type>::setCellIds()")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        V_ = 0.0;
        forAll(cells_, i)
        {
            V_ += mesh_.V()[cells_[i]];
        }
        reduce(V_, sumOp<scalar>());
    }

    Info<< indent << "- selected "
        << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << nl << decrIndent << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeActivatedExplicitSource<Type>::TimeActivatedExplicitSource
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const wordList& fieldNames
)
:
    name_(name),
    mesh_(mesh),
    active_(readBool(dict.lookup("active"))),
    timeStart_(readScalar(dict.lookup("timeStart"))),
    duration_(readScalar(dict.lookup("duration"))),
    volumeMode_(wordToVolumeModeType(dict.lookup("volumeMode"))),
    selectionMode_(wordToSelectionModeType(dict.lookup("selectionMode"))),
    points_(),
    cellSetName_("none"),
    V_(1.0),
    fieldData_(),
    fieldIds_(fieldNames.size(), -1)
{
    setSelection(dict);

    if (fieldNames.size() == 1)
    {
        fieldData_.setSize(1);
        fieldData_[0].first() = fieldNames[0];
        dict.lookup("fieldData") >> fieldData_[0].second();
        fieldIds_[0] = 0;
    }
    else
    {
        setFieldData(dict, fieldNames);
    }

    setCellSet();
}


template<class Type>
void Foam::TimeActivatedExplicitSource<Type>::addToField
(
    DimensionedField<Type, volMesh>& Su,
    const label fieldI
)
{
    const label fid = fieldIds_[fieldI];

    if
    (
        active_
     && (fid >= 0)
     && (mesh_.time().value() >= timeStart_)
     && (mesh_.time().value() <= timeEnd())
    )
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }

        forAll(cells_, i)
        {
            Su[cells_[i]] = fieldData_[fid].second()/V_;
        }
    }
}


// ************************************************************************* //
