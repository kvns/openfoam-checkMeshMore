/*---------------------------------------------------------------------------*\

Class
    Foam::meshChecker

SourceFiles
    meshChecker.C
    meshChecker.H

License
    GPLv3 or later

Copyright
    2013 OpenFOAM Foundation
    2013 Kevin Smith

\*---------------------------------------------------------------------------*/

#include "meshChecker.H"
#include "faceSet.H"
#include "unitConversion.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::meshChecker::meshChecker(const polyMesh& mesh)
:
    _mesh(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshChecker::~meshChecker()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::meshChecker::checkFaceOrthogonality
(
    const bool report,
    labelHashSet* setPtr
) const
{
    Foam::scalar nonOrthThreshold = 70;    // deg
    
    int nBins = 9;
    Foam::scalar binWidth = 10.0;
    List<label> bins(nBins);
    List<Foam::scalar> binVals(nBins+1);
    for (int i = 0;i < nBins; i++)
    {
        bins[i] = 0;
        binVals[i] = ::cos(degToRad(Foam::scalar(i) * binWidth));
    }
    binVals[nBins] = ::cos(degToRad(Foam::scalar(nBins) * binWidth));
    
    // todo, figure out how this debug system works, add back into method
    /*if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceOrthogonality("
            << "const bool, labelHashSet*) const: "
            << "checking mesh non-orthogonality" << endl;
    }*/

    // for all internal faces check that the d dot S product is positive
    const vectorField& centres = _mesh.cellCentres();
    const vectorField& areas = _mesh.faceAreas();

    const labelList& own = _mesh.faceOwner();
    const labelList& nei = _mesh.faceNeighbour();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(nonOrthThreshold));

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll(nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        const vector& s = areas[faceI];

        scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

        for (int binI = 0; binI < nBins; binI++)
        {
            Foam::scalar hi = binVals[binI];
            Foam::scalar low = binVals[binI+1];
            
            if ((low <= dDotS) && (dDotS <= hi))
            {
                bins[binI] += 1;
                break;
            }
        }

        if (dDotS < severeNonorthogonalityThreshold)
        {
            if (dDotS > SMALL)
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                severeNonOrth++;
            }
            else
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                errorNonOrth++;
            }
        }

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }

        sumDDotS += dDotS;
    }
    
    for (int binI = 0; binI < nBins; binI++)
    {
        reduce(bins[binI], sumOp<label>());
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    if ( report)
    {
        label neiSize = nei.size();
        reduce(neiSize, sumOp<label>());
    
        Info<< endl <<"    Mesh non-orthogonality distribution. " <<
               "Binwidth = " << binWidth << " degrees." << endl;
        Info<< "        " <<
             setw(10) << "Angle" << 
             setw(15) << "Number faces" << 
             setw(10) << "Percent faces" << endl;
             
        for (int binI = 0; binI < nBins; binI++)
        {    
            scalar binval(binI * binWidth);
            Info<< "        " << 
                setw(10) << binval << 
                setw(15) << bins[binI] << 
                setw(10) << (bins[binI] / scalar(neiSize))*100 << endl;
        }   

        if (neiSize > 0)
        {
            if (report)
            {
                Info<< "    Mesh non-orthogonality Max: "
                    << radToDeg(::acos(minDDotS))
                    << " average: " << radToDeg(::acos(sumDDotS/neiSize))
                    << endl;
            }
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal "
                << "faces above threshold (" << nonOrthThreshold 
                << " deg): "
                << severeNonOrth << "." << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (report)
        {
            Info<< " ***Non-orthogality errors present. This happens as "
                << "non-orthogonality approaches 90 degrees." << endl;
            Info<< " ***Number of non-orthogonality errors: "
                << errorNonOrth << "." << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Non-orthogonality check OK. No errors detected." 
                << endl;
        }

        return false;
    }
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::meshChecker::operator=(const meshChecker& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::meshChecker::operator=(const Foam::meshChecker&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
