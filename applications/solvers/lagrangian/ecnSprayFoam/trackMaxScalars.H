scalar maxT = max(T).value();
if (TmaxFilePtr.valid())
{
    TmaxFilePtr() << scientific;
    TmaxFilePtr().width(15); TmaxFilePtr() << runTime.timeName() << ",";
    TmaxFilePtr().width(15); TmaxFilePtr() << maxT << endl;
}

label indexOH = -1;
forAll(Y, i)
{
    if (Y[i].name() == "OH")
    {
        indexOH = i;
    }
}
if (indexOH >= 0)
{
    scalar maxOH = max(Y[indexOH]).value();
    if (OHmaxFilePtr.valid())
    {
        OHmaxFilePtr() << scientific;
        OHmaxFilePtr().width(15); OHmaxFilePtr() << runTime.timeName() << ",";
        OHmaxFilePtr().width(15); OHmaxFilePtr() << maxOH << endl;
    }
}
