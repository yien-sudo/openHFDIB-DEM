// Temperature equation
{
    scalar tolEqn = 1e-6;
    label maxEqnIters = 5;
    scalar correctCoeff = 1.0;
    scalar penaltyCoeff = 0.8;
    bool thermalDebug = false;
    const scalar solidMaskThreshold = 0.999;

    if (pimple.dict().found("HFDIB"))
    {
        const dictionary& HFDIBDict = pimple.dict().subDict("HFDIB");
        if (HFDIBDict.found("T"))
        {
            const dictionary& thermalDict = HFDIBDict.subDict("T");
            tolEqn = thermalDict.lookupOrDefault<scalar>("tolEqn", tolEqn);
            maxEqnIters = thermalDict.lookupOrDefault<label>("maxEqnIters", maxEqnIters);
            correctCoeff = thermalDict.lookupOrDefault<scalar>("relaxCoeff", correctCoeff);
            thermalDebug = thermalDict.lookupOrDefault<bool>("thermalDebug", thermalDebug);
            penaltyCoeff = thermalDict.lookupOrDefault<scalar>("penaltyCoeff", penaltyCoeff);
        }
    }

    HFDIBDEM.interpolateScalar(T, Ti, lambda);

    volScalarField Tq
    (
        IOobject
        (
            "Tq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", T.dimensions()/dimTime, 0.0)
    );

    volScalarField TqInterface
    (
        IOobject
        (
            "TqInterface",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", T.dimensions()/dimTime, 0.0)
    );

    const scalarField& lambdaIF = lambda.internalField();

    const scalarField& surfaceIF = surface.internalField();

    const scalar invDeltaT = penaltyCoeff/max(runTime.deltaTValue(), SMALL);

    if (thermalDebug)
    {
        const scalarField& cellVolumes = mesh.V();
        scalar lambdaVolumeLocal = 0.0;
        scalar enforcedVolumeLocal = 0.0;

        forAll(lambdaIF, cellI)
        {
            const scalar cellVol = cellVolumes[cellI];
            const scalar lambdaVal = lambdaIF[cellI];

            lambdaVolumeLocal += lambdaVal*cellVol;
            if (surfaceIF[cellI] > SMALL)
            {
                enforcedVolumeLocal += cellVol;
            }
        }

        scalar lambdaVolume = lambdaVolumeLocal;
        scalar enforcedVolume = enforcedVolumeLocal;

        reduce(lambdaVolume, sumOp<scalar>());
        reduce(enforcedVolume, sumOp<scalar>());

        Info<< "HFDIB thermalDebug: enforcedVol=" << enforcedVolume
            << " lambdaVol=" << lambdaVolume << endl;
    }

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(DT, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);

    for (label iter = 0; iter < maxEqnIters; ++iter)
    {
    const scalarField& TiIF = Ti.internalField();
    const scalarField& TIF = T.internalField();
        const scalarField& cellVolumes = mesh.V();

        scalarField& TqIF = Tq.primitiveFieldRef();
        scalarField& TqInterfaceIF = TqInterface.primitiveFieldRef();

        scalar solidSourceSumLocal = 0.0;
        scalar solidVolumeLocal = 0.0;
        scalar interfaceSourceSumLocal = 0.0;
        scalar interfaceVolumeLocal = 0.0;
        label solidCellCountLocal = 0;
        label interfaceCellCountLocal = 0;

        forAll(TqIF, cellI)
        {
            const scalar surfaceVal = surfaceIF[cellI];
            const bool enforceCell = surfaceVal > SMALL;
            const scalar lambdaVal = lambdaIF[cellI];
            const bool solidCell = lambdaVal >= solidMaskThreshold;
            const scalar interfaceFraction = Foam::max
            (
                scalar(0),
                scalar(1) - lambdaVal
            );
            const bool interfaceCell = enforceCell && (interfaceFraction > SMALL);

            scalar appliedSource = 0.0;
            scalar interfaceContribution = 0.0;

            if (enforceCell)
            {
                const scalar deltaSource = invDeltaT*(TiIF[cellI] - TIF[cellI]);
                appliedSource = surfaceVal*deltaSource;

                if (interfaceCell)
                {
                    interfaceContribution = surfaceVal*interfaceFraction*deltaSource;
                }

                if (thermalDebug)
                {
                    if (solidCell)
                    {
                        solidCellCountLocal++;
                        solidVolumeLocal += lambdaVal*cellVolumes[cellI];
                        solidSourceSumLocal += deltaSource*lambdaVal*cellVolumes[cellI];
                    }
                    else if (interfaceFraction > SMALL)
                    {
                        interfaceCellCountLocal++;
                        interfaceVolumeLocal += interfaceFraction*cellVolumes[cellI];
                        interfaceSourceSumLocal += deltaSource*interfaceFraction*cellVolumes[cellI];
                    }
                }
            }

            TqIF[cellI] = appliedSource;
            TqInterfaceIF[cellI] = interfaceContribution;
        }

        Tq.correctBoundaryConditions();
        TqInterface.correctBoundaryConditions();

        if (thermalDebug)
        {
            scalar solidSourceSum = solidSourceSumLocal;
            scalar solidVolume = solidVolumeLocal;
            scalar interfaceSourceSum = interfaceSourceSumLocal;
            scalar interfaceVolume = interfaceVolumeLocal;
            label solidCellCount = solidCellCountLocal;
            label interfaceCellCount = interfaceCellCountLocal;

            reduce(solidSourceSum, sumOp<scalar>());
            reduce(solidVolume, sumOp<scalar>());
            reduce(interfaceSourceSum, sumOp<scalar>());
            reduce(interfaceVolume, sumOp<scalar>());
            reduce(solidCellCount, sumOp<label>());
            reduce(interfaceCellCount, sumOp<label>());

            Info<< "HFDIB thermalDebug: enforced cells=" << solidCellCount
                << " enforcedVolume=" << solidVolume
                << " interfaceCells=" << interfaceCellCount
                << " interfaceVolume=" << interfaceVolume
                << " solidSourceInt=" << solidSourceSum
                << " interfaceSourceInt=" << interfaceSourceSum
                << endl;
        }

        solve(TEqn == Tq);

    scalar maxError = max(mag(surface*(Ti - T))).value();
        Info<< "HFDIB: max |T - Ti| = " << maxError << endl;

        if (maxError < tolEqn)
        {
            Info<< "HFDIB: T converged to Ti within tolerance " << tolEqn << endl;
            break;
        }

        T += correctCoeff*surface*(Ti - T);
        T.correctBoundaryConditions();
    }

    fvOptions.correct(T);
    T.correctBoundaryConditions();

    HFDIBDEM.updateThermalDiagnostics(T, Tq, TqInterface);
}
