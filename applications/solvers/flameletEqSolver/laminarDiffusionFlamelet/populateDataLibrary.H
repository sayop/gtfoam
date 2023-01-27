{
    forAll(Y, speci)
    {
        //- Given specie's mass fraction in size of entire domain field.
        List<scalar>& data = massFraction[speci];
        //- Re-organize temperature data across all processors
        if (Pstream::parRun())
        {
            //- procData: array of given specie's mass fraction data
            //- in the current processor's decomposed field of domain.
            List<List<scalar>> procData(Pstream::nProcs());
            List<scalar>& myData = procData[Pstream::myProcNo()];
            myData.setSize(Y[speci].size(), 0.0);
            forAll(Y[speci], celli)
            {
                myData[celli] = Y[speci][celli];
            }
            Pstream::gatherList(procData);
 
            if (Pstream::master())
            {
                for (label proci = 0; proci < Pstream::nProcs(); proci++)
                {
                    label dataSize = procData[proci].size();
                    SubList<scalar>
                    (
                        data,
                        dataSize,
                        startIndexList[proci]
                    ) = procData[proci];
                }
            }
            Pstream::scatter(temperature);
            //- Sort data by sortedIndex
            List<scalar> tmp = data;
            forAll(tmp, i)
            {
                tmp[i] = data[sortedIndx[i]];
            }
            data = tmp;
        }
        else
        {
            forAll(Y[speci], celli)
            {
                data[celli] = Y[speci][celli];
            }
        }
    }

    //- Re-organize temperature data across all processors
    if (Pstream::parRun())
    {
        List<List<scalar>> procTemperature(Pstream::nProcs());
        List<scalar>& myTemperature = procTemperature[Pstream::myProcNo()];
        myTemperature.setSize(T.size(), 0.0);
        forAll(T, celli)
        {
            myTemperature[celli] = T[celli];
        }
        Pstream::gatherList(procTemperature);

        if (Pstream::master())
        {
            for (label proci = 0; proci < Pstream::nProcs(); proci++)
            {
                label dataSize = procTemperature[proci].size();
                SubList<scalar>
                (
                    temperature,
                    dataSize,
                    startIndexList[proci]
                ) = procTemperature[proci];
            }
        }
        Pstream::scatter(temperature);
        //- Sort data by sortedIndex
        List<scalar> tmp = temperature;
        forAll(tmp, i)
        {
            tmp[i] = temperature[sortedIndx[i]];
        }
        temperature = tmp;
    }
    else
    {
        forAll(T, celli)
        {
            temperature[celli] = T[celli];
        }
    }

    //- Print out files
    if (Pstream::master())
    {
        forAll(Y, speci)
        {
            OFstream& file = specieDataFilePtrs_[speci];
            file<<runTime.timeName()<<",";
            forAll(massFraction[speci], celli)
            {
                file<<massFraction[speci][celli]<<",";
            }
            file<<endl;
        }

        fileTemperature<<runTime.timeName()<<",";
        forAll(temperature, celli)
        {
            fileTemperature<<temperature[celli]<<",";
        }
        fileTemperature<<endl;
    }
}
