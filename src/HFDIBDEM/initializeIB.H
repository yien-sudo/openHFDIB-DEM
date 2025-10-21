if (mesh_.nGeometricD() < 3)
    {
        const vector validDirs = (mesh_.geometricD() + Vector<label>::one)/2;
        Axis_ -= validDirs;
    }
    Axis_ /= mag(Axis_);

    // read the body operation name from immersedDict_
    if (immersedDict_.found("staticBody"))
    {
        bodyOperation_ = 0;
        InfoH << iB_Info << bodyName_ << " is a static body." << endl;
    }
    else if (immersedDict_.found("prescribedTransBody"))
    {
        bodyOperation_ = 1;

        Vel_   = immersedDict_.subDict("prescribedTransBody").lookup("velocity");

        InfoH << iB_Info << bodyName_ << " is a freely rotating body with prescribed linear velocity." << endl;
    }
    else if (immersedDict_.found("prescribedRotBody"))
    {
        bodyOperation_ = 2;

        Axis_ = immersedDict_.subDict("prescribedRotBody").lookup("axis");
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedRotBody").lookup("omega")
                );

        InfoH << iB_Info << bodyName_ << " is a freely moving body with prescribed rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransRotBody"))
    {
        bodyOperation_ = 3;

        Vel_   = immersedDict_.subDict("prescribedTransRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransRotBody").lookup("axis");
        omega_  = readScalar(
                    immersedDict_.subDict("prescribedTransRotBody").lookup("omega")
                );

        InfoH << iB_Info << bodyName_ << " has prescribed both movement and rotation." << endl;
    }
    else if (immersedDict_.found("prescribedTransFixedAxisRotBody"))
    {
        bodyOperation_ = 4;

        Vel_   = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("velocity");
        Axis_ = immersedDict_.subDict("prescribedTransFixedAxisRotBody").lookup("axis");

        InfoH << iB_Info << bodyName_ << " has prescribed movement and axis of rotation." << endl;
    }
    else if (immersedDict_.found("fullyCoupledBody"))
    {
        bodyOperation_ = 5;

        if (immersedDict_.isDict("fullyCoupledBody"))
        {
            if(immersedDict_.subDict("fullyCoupledBody").found("velocity"))
            {
                Vel_ = immersedDict_.subDict("fullyCoupledBody").lookup("velocity");
            }
            if(immersedDict_.subDict("fullyCoupledBody").found("axis"))
            {            
                Axis_ = immersedDict_.subDict("fullyCoupledBody").lookup("axis");
            }
            if(immersedDict_.subDict("fullyCoupledBody").found("omega"))
            {
                omega_  = readScalar(
                    immersedDict_.subDict("fullyCoupledBody").lookup("omega")
                );
            }
        }

        InfoH << iB_Info << bodyName_ << " is fully coupled with fluid phase." << endl;
    }
    else
    {
        InfoH << iB_Info << "No body operation was found for " << bodyName_ << endl
          << "Assuming static body.";
    }

    // check if the immersedDict_ contains switch for updateTorque_
    if (immersedDict_.found("updateTorque"))
    {
        updateTorque_ = readBool(immersedDict_.lookup("updateTorque"));
        InfoH << iB_Info << "Found updateTorque: " << updateTorque_ << endl;
    }
    else
    {
        InfoH << iB_Info << "Did not find updateTorque, using updateTorque: " << updateTorque_ << endl;
    }

    // do I want the IB to start as in sync with the flow?
    if (immersedDict_.found("startSynced"))
    {
        startSynced_ = readBool(immersedDict_.lookup("startSynced"));
        if (startSynced_)
        {
            InfoH << iB_Info << "Will try to sync the body with the flow upon creation" << endl;
        }
        else
        {
            InfoH << iB_Info << "The body will be created as static" << endl;
        }
    }
    else
    {
        InfoH << iB_Info << "startSynced was not specified, using startSynced: " << startSynced_ << endl;
    }

    // body auxiliary settings
    if (immersedDict_.found("sdBasedLambda"))
    {
        sdBasedLambda_ = readBool(immersedDict_.lookup("sdBasedLambda"));
    }
    if (immersedDict_.found("interfaceSpan"))
    {
        intSpan_ = readScalar(immersedDict_.lookup("interfaceSpan"));
    }
    if (immersedDict_.found("refineBuffers"))
    {
        refineBuffers_ = readLabel(immersedDict_.lookup("refineBuffers"));
    }
    if (immersedDict_.found("timesToSetStatic"))
    {
        timesToSetStatic_ = readLabel(immersedDict_.lookup("timesToSetStatic"));
    }

    // set up body geometry
    geomModel_->setIntSpan(intSpan_);
    geomModel_->setSdBasedLambda(sdBasedLambda_);
    geomModel_->setRhoS(immersedDict_.lookup("rho"));

    string material = word(immersedDict_.lookup("material"));

    if(ibIntp.valid())
    {
        intpInfo_.set(ibIntp().createIntpInfo(
            mesh_,
            geomModel_
        ));
    }

    ibContactClass_.set(new ibContactClass(
        geomModel_,
        material
    ));

    ibCntVars_.set(new ibContactVars
    (
        bodyId_,
        Vel_,
        omega_,
        Axis_,
        geomModel_->getM0(),
        geomModel_->getM(),
        geomModel_->getRhoS()
    ));

    wallCntInfo_.set(new wallContactInfo
    (
        ibContactClass_(),
        ibCntVars_()
    ));

    if (immersedDict_.found("T"))
    {
        hasThermalBC_ = true;

        const dictionary& thermalDict = immersedDict_.subDict("T");

        scalar boundaryValue = 0.0;
        if (thermalDict.found("value"))
        {
            boundaryValue = readScalar(thermalDict.lookup("value"));
        }

        scalar innerValue = boundaryValue;
        if (thermalDict.found("innerValue"))
        {
            innerValue = readScalar(thermalDict.lookup("innerValue"));
        }

        surfaceTemperature_ = boundaryValue;
        particleTemperature_ = innerValue;
        if (thermalDict.found("bulkTemperature"))
        {
            bulkTemperature_ = readScalar(thermalDict.lookup("bulkTemperature"));
            bulkTemperatureFixed_ = true;
        }
        else
        {
            bulkTemperature_ = boundaryValue;
            bulkTemperatureFixed_ = false;
        }

        thermalDynamic_ = thermalDict.lookupOrDefault<bool>("dynamic", false);
        computeNu_ = thermalDict.lookupOrDefault<bool>("computeNu", false);
    thermalDebug_ = thermalDict.lookupOrDefault<bool>("thermalDebug", false);

        cpParticle_ = thermalDict.lookupOrDefault<scalar>("cpParticle", 0.0);
        cpFluid_ = thermalDict.lookupOrDefault<scalar>("cpFluid", 0.0);
        kFluid_ = thermalDict.lookupOrDefault<scalar>("lambdaFluid", 0.0);

        referenceDiameter_ = thermalDict.lookupOrDefault<scalar>
        (
            "characteristicLength",
            getDC()
        );

        referenceArea_ = thermalDict.lookupOrDefault<scalar>(
            "referenceArea",
            0.0
        );

        minDeltaT_ = Foam::max
        (
            thermalDict.lookupOrDefault<scalar>("minDeltaT", minDeltaT_),
            SMALL
        );

        solidFractionThreshold_ = Foam::max
        (
            Foam::min
            (
                thermalDict.lookupOrDefault<scalar>
                (
                    "solidFractionThreshold",
                    solidFractionThreshold_
                ),
                1.0
            ),
            0.0
        );

        if (kFluid_ <= SMALL && cpFluid_ > SMALL)
        {
            kFluid_ = rhoF_.value()*cpFluid_*DTF_.value();
        }
        else if (cpFluid_ <= SMALL && kFluid_ > SMALL)
        {
            cpFluid_ = kFluid_/(rhoF_.value()*DTF_.value());
        }

        thermalConfigAvailable_ = (cpFluid_ > SMALL);

        if (!thermalConfigAvailable_)
        {
            if (thermalDynamic_ || computeNu_)
            {
                InfoH << iB_Info
                      << "Thermal properties for body " << bodyName_
                      << " are incomplete (cpFluid missing)."
                      << " Thermal forcing diagnostics disabled." << endl;
            }
            thermalDynamic_ = false;
            computeNu_ = false;
        }

        if (thermalDynamic_ && cpParticle_ <= SMALL)
        {
            InfoH << iB_Info
                  << "Particle heat capacity cpParticle missing for body "
                  << bodyName_ << ". Disabling thermal evolution." << endl;
            thermalDynamic_ = false;
        }

        if (computeNu_ && kFluid_ <= SMALL)
        {
            InfoH << iB_Info
                  << "Thermal conductivity kFluid missing for body "
                  << bodyName_ << ". Nusselt number output disabled." << endl;
            computeNu_ = false;
        }

        if (referenceDiameter_ <= SMALL)
        {
            const scalar rhoSolid = geomModel_->getRhoS().value();
            const scalar massSolid = geomModel_->getM0();
            if (rhoSolid > SMALL && massSolid > SMALL)
            {
                const scalar volume = massSolid/rhoSolid;
                if (volume > SMALL)
                {
                    referenceDiameter_ = std::cbrt(
                        6.0*volume/(Foam::constant::mathematical::pi)
                    );
                }
            }
        }

        if (referenceArea_ <= SMALL && referenceDiameter_ > SMALL)
        {
            referenceArea_ = Foam::constant::mathematical::pi
                *sqr(referenceDiameter_);
        }
    }
