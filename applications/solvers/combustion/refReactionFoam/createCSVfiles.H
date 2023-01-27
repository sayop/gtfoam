autoPtr<OFstream> TmaxFilePtr;
autoPtr<OFstream> OHmaxFilePtr;
if (Pstream::master())
{
    TmaxFilePtr.reset
    (
        new OFstream("Tmax.csv")
    );
    TmaxFilePtr().width(15); TmaxFilePtr() << "Time [sec]" << ",";
    TmaxFilePtr().width(15); TmaxFilePtr() << "max(Temp) [K]" << endl;

    OHmaxFilePtr.reset
    (
        new OFstream("OHmax.csv")
    );
    OHmaxFilePtr().width(15); OHmaxFilePtr() << "Time [sec]" << ",";
    OHmaxFilePtr().width(15); OHmaxFilePtr() << "max(OH) [-]" << endl;

};


