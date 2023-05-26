import numpy as np
from scipy.integrate import solve_ivp 
from .axis import Axis
from mpy.specMagneticField import FieldLine, specField
from mpy.specMagneticField import readJacobian
from typing import Tuple


def findBifurcation(firstAxis: Axis, secondAxis: Axis, bField: specField, jacobianData: str, niter: int=10, plotDebug: bool=False, iterLine: int=6) -> Tuple[float]:
    
    if jacobianData is None:
        base_Jacobian = bField.getB()
        base_sArr = bField.sArr
        base_thetaArr = bField.thetaArr
        base_zetaArr = bField.zetaArr
    else:
        base_sArr, base_thetaArr, base_zetaArr, base_Jacobian = readJacobian(jacobianData)

    firstR = float(firstAxis.getRZ(np.array([0]))[0])
    secondR = float(secondAxis.getRZ(np.array([0]))[0])
    print("R = " + "{:.2e}".format(firstR) + ", " + "{:.2e}".format(secondR))
    assert firstR < secondR
    leftS, rightS = firstAxis.initPoint[0], secondAxis.initPoint[0]
    # leftS, rightS = -1, 1
    # leftState, rightState = True, False

    if plotDebug:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()

    lineArr = list()
    for _i in range(niter):
        midS = (leftS + rightS) / 2
        lineArr.append(traceLine(np.array([midS, 0]), bField, base_sArr, base_thetaArr, base_zetaArr, base_Jacobian, iterLine))
        rArr = list()
        zArr = list()
        for i in range(len(lineArr[-1].rArr)):
            if i % lineArr[-1].nZeta == 0:
                rArr.append(lineArr[-1].rArr[i])
                zArr.append(lineArr[-1].zArr[i])
        if plotDebug:
            dots = ax.scatter(rArr, zArr, s=2.0)
        midR = max(rArr)
        print("midR = " + "{:.2e}".format(midR))
        assert midR > firstR
        if midR < secondR:
            leftS = midS
        else:
            rightS = midS

    if plotDebug:
        plt.axis("equal")
    
    return leftS, rightS


def traceLine(initPoint: np.ndarray, bField: specField, 
    base_sArr: np.ndarray, base_thetaArr: np.ndarray, base_zetaArr: np.ndarray, base_Jacobian: np.ndarray, 
    iterLine: int, nstep: int=4) -> FieldLine:
    
    import pyoculus
    pyoculusField = pyoculus.problems.SPECBfield(bField.specData, bField.lvol+1)
    def getB(zeta, s_theta):
        field = pyoculusField.B_many(s_theta[0], s_theta[1], zeta) / bField.interpValue(base_Jacobian, s_theta[0], s_theta[1], zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        bSupS = field[0, 0]
        bSupTheta = field[0, 1]
        bSupZeta = field[0, 2]
        return [bSupS/bSupZeta, bSupTheta/bSupZeta]
    sValue, thetaValue = initPoint
    s_theta = [sValue, thetaValue]
    zetaStart = 0
    dZeta = 2 * np.pi / bField.nfp / nstep
    sArr = [sValue]
    thetaArr = [thetaValue]
    zetaArr = [0]
    niter = bField.nfp * iterLine
    for j in range(niter):
        for k in range(nstep):
            sol = solve_ivp(getB, (zetaStart,zetaStart+dZeta), s_theta, method="LSODA", rtol=1e-9)
            sArr.append(sol.y[0,-1])
            thetaArr.append(sol.y[1,-1])
            zetaArr.append(zetaStart+dZeta)
            s_theta = [sArr[-1], thetaArr[-1]]
            zetaStart = zetaArr[-1]
    line = FieldLine.getLine_tracing(bField, nstep, np.array(sArr), np.array(thetaArr), np.array(zetaArr))
    return line



if __name__ == "__main__":
    pass
