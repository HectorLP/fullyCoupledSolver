#include "fvCFD.H"

scalar computeR(const fvMesh& mesh, volScalarField& r, dimensionedVector x0);

void computeU(const fvMesh& mesh, volVectorField& U, word pName = "p");
