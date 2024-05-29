#include "customLibrary.H"

scalar computeR(const fvMesh& mesh, volScalarField& r, dimensionedVector x0)
{
    r = mag(mesh.C() - x0);
    return returnReduce(max(r).value(), maxOp<scalar>());
}

void computeU(const fvMesh& mesh, volVectorField& U, word pName)
{
    const volScalarField& pField = mesh.lookupObject<volScalarField>(pName);

    U = fvc::grad(pField)*dimensionedScalar("tmp", dimTime, 1.);
}
