/*---------------------------------------------------------------------------*\

Class
    Foam::meshChecker

Description
    Class wraps modified mesh quality checking to include additional information
    about a given mesh.

SourceFiles
    meshChecker.C
    meshChecker.H

License: 
    GPLv3 or later

Copyright:
    2013 Kevin Smith

\*---------------------------------------------------------------------------*/

#ifndef meshChecker_H
#define meshChecker_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshChecker Declaration
\*---------------------------------------------------------------------------*/

class meshChecker

{
    // Private data
    const polyMesh& _mesh;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshChecker(const meshChecker&);

        //- Disallow default bitwise assignment
        void operator=(const meshChecker&);

        

public:


    // Constructors
        
        meshChecker(const polyMesh&);

    //- Destructor
    ~meshChecker();


    // Member Functions
    bool checkFaceOrthogonality(
            const bool report, 
            labelHashSet* setPtr) const;
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
