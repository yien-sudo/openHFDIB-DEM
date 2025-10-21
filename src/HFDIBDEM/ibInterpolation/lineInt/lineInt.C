/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "lineInt.H"

#define ORDER 2

using namespace Foam;

namespace
{
template<class Type, class Setter>
void interpolateIntPointsParallel
(
    const fvMesh& mesh,
    interpolation<Type>& interpolator,
    List<List<intPoint>>& intPoints,
    Setter&& setter
)
{
    const label myProc = Pstream::myProcNo();
    const label nCells = mesh.nCells();

    List<DynamicList<Tuple2<label, label>>> owners(Pstream::nProcs());
    List<DynamicList<point>> pointRequests(Pstream::nProcs());
    List<DynamicList<label>> cellRequests(Pstream::nProcs());

    forAll(intPoints, ibCellI)
    {
        List<intPoint>& intList = intPoints[ibCellI];

        forAll(intList, iPoint)
        {
            intPoint& curPoint = intList[iPoint];

            setter(curPoint, Type());

            if (curPoint.iProc_ == myProc)
            {
                if (curPoint.iCell_ < 0 || curPoint.iCell_ >= nCells)
                {
                    curPoint.iProc_ = -1;
                    continue;
                }

                setter
                (
                    curPoint,
                    interpolator.interpolate(curPoint.iPoint_, curPoint.iCell_)
                );
            }
            else if (curPoint.iProc_ != -1)
            {
                pointRequests[curPoint.iProc_].append(curPoint.iPoint_);
                cellRequests[curPoint.iProc_].append(curPoint.iCell_);
                owners[curPoint.iProc_].append
                (
                    Tuple2<label, label>(ibCellI, iPoint)
                );
            }
        }
    }

    PstreamBuffers pointBufs(Pstream::commsTypes::nonBlocking);
    PstreamBuffers cellBufs(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci == myProc)
        {
            continue;
        }

        UOPstream sendPoints(proci, pointBufs);
        UOPstream sendCells(proci, cellBufs);

        sendPoints << pointRequests[proci];
        sendCells << cellRequests[proci];
    }

    pointBufs.finishedSends();
    cellBufs.finishedSends();

    List<DynamicList<point>> pointRecv(Pstream::nProcs());
    List<DynamicList<label>> cellRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci == myProc)
        {
            continue;
        }

        UIPstream recvPoints(proci, pointBufs);
        UIPstream recvCells(proci, cellBufs);

    pointRecv[proci] = DynamicList<point>(recvPoints);
    cellRecv[proci] = DynamicList<label>(recvCells);
    }

    pointBufs.clear();
    cellBufs.clear();

    List<DynamicList<Type>> valueReturn(Pstream::nProcs());

    forAll(pointRecv, proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        forAll(pointRecv[proci], idx)
        {
            const label cellI = cellRecv[proci][idx];

            if (cellI < 0 || cellI >= nCells)
            {
                valueReturn[proci].append(Type());
                continue;
            }

            valueReturn[proci].append
            (
                interpolator.interpolate(pointRecv[proci][idx], cellI)
            );
        }
    }

    PstreamBuffers valueBufs(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci == myProc)
        {
            continue;
        }

        UOPstream sendVals(proci, valueBufs);
        sendVals << valueReturn[proci];
    }

    valueBufs.finishedSends();

    List<DynamicList<Type>> valueRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci == myProc)
        {
            continue;
        }

        UIPstream recvVals(proci, valueBufs);
        valueRecv[proci] = DynamicList<Type>(recvVals);
    }

    valueBufs.clear();

    forAll(valueRecv, proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        const DynamicList<Tuple2<label, label>>& ownerList = owners[proci];

        forAll(valueRecv[proci], idx)
        {
            if (idx >= ownerList.size())
            {
                break;
            }

            const Tuple2<label, label>& owner = ownerList[idx];

            if
            (
                owner.first() < 0
             || owner.first() >= intPoints.size()
            )
            {
                continue;
            }

            List<intPoint>& slot = intPoints[owner.first()];

            if
            (
                owner.second() < 0
             || owner.second() >= slot.size()
            )
            {
                continue;
            }

            setter(slot[owner.second()], valueRecv[proci][idx]);
        }
    }
}

template<class Setter>
void gatherCellFieldParallel
(
    const fvMesh& mesh,
    const scalarField& cellField,
    List<List<intPoint>>& intPoints,
    Setter&& setter
)
{
    const label myProc = Pstream::myProcNo();
    const label nCells = mesh.nCells();

    List<DynamicList<Tuple2<label, label>>> owners(Pstream::nProcs());
    List<DynamicList<label>> cellRequests(Pstream::nProcs());

    forAll(intPoints, ibCellI)
    {
        List<intPoint>& intList = intPoints[ibCellI];

        forAll(intList, iPoint)
        {
            intPoint& curPoint = intList[iPoint];

            if (curPoint.iProc_ == myProc)
            {
                if (curPoint.iCell_ >= 0 && curPoint.iCell_ < nCells)
                {
                    setter(curPoint, cellField[curPoint.iCell_]);
                }
            }
            else if (curPoint.iProc_ != -1)
            {
                cellRequests[curPoint.iProc_].append(curPoint.iCell_);
                owners[curPoint.iProc_].append
                (
                    Tuple2<label, label>(ibCellI, iPoint)
                );
            }
        }
    }

    PstreamBuffers cellBufs(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        UOPstream sendCells(proci, cellBufs);
        sendCells << cellRequests[proci];
    }

    cellBufs.finishedSends();

    List<DynamicList<label>> cellRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        UIPstream recvCells(proci, cellBufs);
        cellRecv[proci] = DynamicList<label>(recvCells);
    }

    cellBufs.clear();

    List<DynamicList<scalar>> valueReturn(Pstream::nProcs());

    forAll(cellRecv, proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        const DynamicList<label>& requests = cellRecv[proci];
        DynamicList<scalar>& responses = valueReturn[proci];

        forAll(requests, idx)
        {
            const label cellI = requests[idx];

            scalar response = 1.0;
            if (cellI >= 0 && cellI < nCells)
            {
                response = cellField[cellI];
            }

            responses.append(response);
        }
    }

    PstreamBuffers valueBufs(Pstream::commsTypes::nonBlocking);

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        UOPstream sendVals(proci, valueBufs);
        sendVals << valueReturn[proci];
    }

    valueBufs.finishedSends();

    List<DynamicList<scalar>> valueRecv(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        UIPstream recvVals(proci, valueBufs);
        valueRecv[proci] = DynamicList<scalar>(recvVals);
    }

    valueBufs.clear();

    forAll(valueRecv, proci)
    {
        if (proci == myProc)
        {
            continue;
        }

        const DynamicList<Tuple2<label, label>>& ownerList = owners[proci];
        const DynamicList<scalar>& values = valueRecv[proci];

        forAll(values, idx)
        {
            if (idx >= ownerList.size())
            {
                break;
            }

            const Tuple2<label, label>& owner = ownerList[idx];

            if
            (
                owner.first() < 0
             || owner.first() >= intPoints.size()
            )
            {
                continue;
            }

            List<intPoint>& slot = intPoints[owner.first()];

            if
            (
                owner.second() < 0
             || owner.second() >= slot.size()
            )
            {
                continue;
            }

            setter(slot[owner.second()], values[idx]);
        }
    }
}
}

//---------------------------------------------------------------------------//
lineInt::lineInt(dictionary& interpDict):
interpDict_(interpDict),
fluidFractionThreshold_
(
    Foam::max
    (
        scalar(0),
        Foam::min
        (
            scalar(1),
            interpDict.lookupOrDefault<scalar>("minFluidFraction", 0.5)
        )
    )
)
{}
lineInt::~lineInt()
{}
//---------------------------------------------------------------------------//
void lineInt::ibInterpolate
(
    interpolationInfo& intpInfo,
    volVectorField& Ui,
    vectorField ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    lineIntInfo& lsInfo
        = dynamic_cast<lineIntInfo&>(intpInfo);

    correctVelocity
    (
        lsInfo,
        Ui,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void lineInt::ibInterpolateScalar
(
    interpolationInfo& intpInfo,
    volScalarField& Si,
    const scalarField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    lineIntInfo& lsInfo = dynamic_cast<lineIntInfo&>(intpInfo);

    correctScalar
    (
        lsInfo,
        Si,
        ibPointsVal,
        mesh
    );
}
//---------------------------------------------------------------------------//
void lineInt::correctVelocity
(
    lineIntInfo& intpInfo,
    volVectorField& Ui,
    vectorField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& cSurfCells = intpInfo.getSurfCells();

    List<point>& ibPoints = intpInfo.getIbPoints();
    List<List<intPoint>>& intPoints = intpInfo.getIntPoints();

    getCurVelocity(intPoints, mesh);
    updateFluidFractions(intPoints, mesh);
    List<label> intOrder = getIntOrder(intPoints);

    forAll(intPoints, ibp)
    {
        label cellI = cSurfCells[ibp];

        switch(intOrder[ibp])
        {
            case 0:
            {
                Ui[cellI] = ibPointsVal[ibp];
                break;
            }

            case 1:
            {
                vector VP1 = intPoints[ibp][0].iVel_ - ibPointsVal[ibp];

                // distance between interpolation points
                scalar deltaR = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);

                // cell center to surface distance
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector linCoeff = VP1/(deltaR+SMALL);

                Ui[cellI] = linCoeff*ds + ibPointsVal[ibp];
                break;
            }

            case 2:
            {
                vector VP1 =  intPoints[ibp][0].iVel_ - ibPointsVal[ibp];
                vector VP2 =  intPoints[ibp][1].iVel_ - ibPointsVal[ibp];

                // distance between interpolation points
                scalar deltaR1 = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);
                scalar deltaR2 = mag(intPoints[ibp][1].iPoint_
                        - intPoints[ibp][0].iPoint_);

                // cell center to surface distance
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                vector quadCoeff = (VP2 - VP1)*deltaR1 - VP1*deltaR2;
                quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                vector linCoeff  = (VP1-VP2)*Foam::pow(deltaR1,2.0);
                linCoeff        += 2.0*VP1*deltaR1*deltaR2;
                linCoeff        += VP1*Foam::pow(deltaR2,2.0);
                linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2)+SMALL);

                Ui[cellI] = quadCoeff*ds*ds + linCoeff*ds + ibPointsVal[ibp];
                break;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void lineInt::correctScalar
(
    lineIntInfo& intpInfo,
    volScalarField& Si,
    const scalarField& ibPointsVal,
    const Foam::fvMesh& mesh
)
{
    const DynamicLabelList& cSurfCells = intpInfo.getSurfCells();

    List<point>& ibPoints = intpInfo.getIbPoints();
    List<List<intPoint>>& intPoints = intpInfo.getIntPoints();

    getCurScalar(intPoints, mesh);
    updateFluidFractions(intPoints, mesh);
    List<label> intOrder = getIntOrder(intPoints);

    if (ibPointsVal.empty())
    {
        return;
    }

    forAll(intPoints, ibp)
    {
        label cellI = cSurfCells[ibp];

        switch(intOrder[ibp])
        {
            case 0:
            {
                Si[cellI] = ibPointsVal[ibp];
                break;
            }

            case 1:
            {
                scalar SP1 = intPoints[ibp][0].iScalar_ - ibPointsVal[ibp];

                scalar deltaR = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                scalar linCoeff = SP1/(deltaR + SMALL);

                Si[cellI] = linCoeff*ds + ibPointsVal[ibp];
                break;
            }

            case 2:
            {
                scalar SP1 = intPoints[ibp][0].iScalar_ - ibPointsVal[ibp];
                scalar SP2 = intPoints[ibp][1].iScalar_ - ibPointsVal[ibp];

                scalar deltaR1 = mag(intPoints[ibp][0].iPoint_ - ibPoints[ibp]);
                scalar deltaR2 = mag(intPoints[ibp][1].iPoint_ - intPoints[ibp][0].iPoint_);
                scalar ds = mag(mesh.C()[cellI] - ibPoints[ibp]);

                scalar quadCoeff = (SP2 - SP1)*deltaR1 - SP1*deltaR2;
                quadCoeff       /= (deltaR1*deltaR2*(deltaR1 + deltaR2) + SMALL);

                scalar linCoeff  = (SP1 - SP2)*Foam::pow(deltaR1, 2.0);
                linCoeff        += 2.0*SP1*deltaR1*deltaR2;
                linCoeff        += SP1*Foam::pow(deltaR2, 2.0);
                linCoeff        /= (deltaR1*deltaR2*(deltaR1 + deltaR2) + SMALL);

                Si[cellI] = quadCoeff*ds*ds + linCoeff*ds + ibPointsVal[ibp];
                break;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void lineInt::getCurVelocity
(
    List<List<intPoint>>& intPoints,
    const fvMesh& mesh
)
{
    if (Pstream::master())
    {
        Info<< "lineInt::getCurVelocity start" << endl;
    }

    if (!interpV_.valid())
    {
        label pending = 0;

        forAll(intPoints, ibCellI)
        {
            forAll(intPoints[ibCellI], iPoint)
            {
                const intPoint& curPoint = intPoints[ibCellI][iPoint];

                if
                (
                    curPoint.iProc_ != -1
                 && curPoint.iProc_ != Pstream::myProcNo()
                )
                {
                    ++pending;
                }
            }
        }

        if (pending)
        {
            Info<< "lineInt::getCurVelocity proc " << Pstream::myProcNo()
                << " missing velocity interpolator with " << pending
                << " remote requests" << endl;
        }

        return;
    }

    interpolateIntPointsParallel<vector>
    (
        mesh,
        *interpV_,
        intPoints,
        [](intPoint& target, const vector& value)
        {
            target.iVel_ = value;
        }
    );

    if (Pstream::master())
    {
        Info<< "lineInt::getCurVelocity end" << endl;
    }
}
//---------------------------------------------------------------------------//
void lineInt::getCurScalar
(
    List<List<intPoint>>& intPoints,
    const fvMesh& mesh
)
{
    if (Pstream::master())
    {
        Info<< "lineInt::getCurScalar start" << endl;
    }

    if (!interpS_.valid())
    {
        label pending = 0;

        forAll(intPoints, ibCellI)
        {
            forAll(intPoints[ibCellI], iPoint)
            {
                const intPoint& curPoint = intPoints[ibCellI][iPoint];

                if
                (
                    curPoint.iProc_ != -1
                 && curPoint.iProc_ != Pstream::myProcNo()
                )
                {
                    ++pending;
                }
            }
        }

        if (pending)
        {
            Info<< "lineInt::getCurScalar proc " << Pstream::myProcNo()
                << " missing scalar interpolator with " << pending
                << " remote requests" << endl;
        }
        return;
    }

    interpolateIntPointsParallel<scalar>
    (
        mesh,
        *interpS_,
        intPoints,
        [](intPoint& target, const scalar value)
        {
            target.iScalar_ = value;
        }
    );

    if (Pstream::master())
    {
        Info<< "lineInt::getCurScalar end" << endl;
    }
}
//---------------------------------------------------------------------------//
void lineInt::updateFluidFractions
(
    List<List<intPoint>>& intPoints,
    const fvMesh& mesh
)
{
    forAll(intPoints, ibp)
    {
        List<intPoint>& pointList = intPoints[ibp];
        forAll(pointList, ip)
        {
            pointList[ip].fluidFraction_ = 1.0;
        }
    }

    if (!mesh.foundObject<volScalarField>("lambda"))
    {
        return;
    }

    const volScalarField& lambdaField =
        mesh.lookupObject<volScalarField>("lambda");
    const scalarField& lambdaIF = lambdaField.internalField();

    forAll(intPoints, ibp)
    {
        List<intPoint>& pointList = intPoints[ibp];
        forAll(pointList, ip)
        {
            pointList[ip].fluidFraction_ = 0.0;
        }
    }

    gatherCellFieldParallel
    (
        mesh,
        lambdaIF,
        intPoints,
        [](intPoint& target, const scalar lambdaVal)
        {
            const scalar clipped = Foam::min
            (
                scalar(1),
                Foam::max(scalar(0), lambdaVal)
            );
            target.fluidFraction_ = scalar(1) - clipped;
        }
    );
}
//---------------------------------------------------------------------------//
List<label> lineInt::getIntOrder
(
    const List<List<intPoint>>& intPoints
) const
{
    List<label> intOrderToRtn(intPoints.size(), 0);

    forAll(intPoints, ibp)
    {
        const List<intPoint>& pointList = intPoints[ibp];

        forAll(pointList, ip)
        {
            const intPoint& curPoint = pointList[ip];

            if
            (
                curPoint.iProc_ != -1
             && curPoint.fluidFraction_ > fluidFractionThreshold_
            )
            {
                ++intOrderToRtn[ibp];
            }
        }
    }

    return intOrderToRtn;
}
//---------------------------------------------------------------------------//
