/*---------------------------------------------------------------------------*\

Application
    checkMeshMore

Description
    Checks validity of a mesh. 
    This application is based on the original checkMesh with additional 
    output information about the mesh and the parameters by which the mesh 
    quality is quantified

License: 
    GPLv3 or later

Copyright:
    2013 OpenFOAM Foundation
    2013 Kevin Smith

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

#include "polyMesh.H"
#include "globalMeshData.H"

#include "printMeshStats.H"
#include "checkTopology.H"
#include "checkGeometry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    Foam::argList::addOption
    (
        "region",
        "name",
        "specify alternative mesh region"
    );

    argList::addBoolOption
    (
        "noTopology",
        "skip checking the mesh topology"
    );
    argList::addBoolOption
    (
        "allGeometry",
        "include bounding box checks"
    );
    argList::addBoolOption
    (
        "allTopology",
        "include extra topology checks"
    );

    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }


    Foam::Info<< "Create time\n" << Foam::endl;

    Foam::Time runTime(Foam::Time::controlDictName, args);


    instantList timeDirs = timeSelector::select0(runTime, args);

    Foam::word regionName;

    if (args.optionReadIfPresent("region", regionName))
    {
        Foam::Info
            << "Create polyMesh " << regionName << " for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }
    else
    {
        regionName = Foam::polyMesh::defaultRegion;
        Foam::Info
            << "Create polyMesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }

    Foam::polyMesh mesh
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    const bool noTopology  = args.optionFound("noTopology");
    const bool allGeometry = args.optionFound("allGeometry");
    const bool allTopology = args.optionFound("allTopology");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        polyMesh::readUpdateState state = mesh.readUpdate();

        if
        (
            !timeI
         || state == polyMesh::TOPO_CHANGE
         || state == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Clear mesh before checking
            mesh.clearOut();

            // Reconstruct globalMeshData
            mesh.globalData();

            printMeshStats(mesh, allTopology);

            label noFailedChecks = 0;

            if (!noTopology)
            {
                noFailedChecks += checkTopology(mesh, allTopology, allGeometry);
            }

            noFailedChecks += checkGeometry(mesh, allGeometry);

            // Note: no reduction in noFailedChecks necessary since is
            //       counter of checks, not counter of failed cells,faces etc.


            if (noFailedChecks == 0)
            {
                Info<< "\nMesh OK.\n" << endl;
            }
            else
            {
                Info<< "\nFailed " << noFailedChecks << " mesh checks.\n"
                    << endl;
            }
        }
        else if (state == polyMesh::POINTS_MOVED)
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            label nFailedChecks = checkGeometry(mesh, allGeometry);

            if (nFailedChecks)
            {
                Info<< "\nFailed " << nFailedChecks << " mesh checks.\n"
                    << endl;
            }
            else
            {
                Info<< "\nMesh OK.\n" << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
