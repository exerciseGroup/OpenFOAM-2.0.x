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

#include "chemkinReader.H"
#include "argList.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("CHEMKINIIIFile");
    argList::addOption("thermo", "fileName");
    argList args(argc, argv);

    fileName thermoFileName = fileName::null;
    args.optionReadIfPresent("thermo", thermoFileName);

    chemkinReader ck(args[1], thermoFileName);

    //Info<< ck.isotopeAtomicWts() << nl
    //    << ck.specieNames() << nl
    //    << ck.speciePhase() << nl
    //    << ck.specieThermo() << nl
    //    << ck.reactions() << endl;

    const SLPtrList<gasReaction>& reactions = ck.reactions();

    {
        OFstream reactionStream("reactions");
        reactionStream<< reactions << endl;
    }

    {
        IFstream reactionStream("reactions");

        label nReactions(readLabel(reactionStream));
        reactionStream.readBeginList(args.executable().c_str());

        PtrList<gasReaction> testReactions(nReactions);

        forAll(testReactions, i)
        {
            testReactions.set
            (
                i,
                gasReaction::New
                (
                    ck.species(),
                    ck.speciesThermo(),
                    reactionStream
                )
            );
        }

        reactionStream.readEndList(args.executable().c_str());

        Info<< testReactions << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
