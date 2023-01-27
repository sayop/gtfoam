combustion->updateSpecie(Z, Zvar, Y);
//label inertIndex = -1;
//volScalarField Yt(0.0*Y[0]);

//forAll(Y, i)
//{
//    combustion->updateSpecie(Y[i]);
//    if (Y[i].name() != inertSpecie)
//    {
//        volScalarField& Yi = Y[i];
//
//        //update Y from flamelet libraries
//        
//        Yi.max(0.0);
//        Yt += Yi;
//    }
//    else
//    {
//        inertIndex = i;
//    }
//}
//Y[inertIndex] = scalar(1) - Yt;
//Y[inertIndex].max(0.0);
