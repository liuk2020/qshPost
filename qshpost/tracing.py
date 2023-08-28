import numpy as np
from scipy.integrate import solve_ivp 
from mpy.specMagneticField import SPECField, FieldLine
from mpy.misc import print_progress
from typing import List


def traceLine(
    bField: SPECField, 
    s0: np.ndarray, theta0: np.ndarray, zeta0: np.ndarray, 
    oneLength: float, 
    niter: int=128, nstep: int=32, 
    sResolution: int=128, thetaResolution: int=128, zetaResolution: int=128, 
    **kwargs
) -> List[FieldLine]:
    r"""
    Working in SPEC coordintes (s, \theta, \zeta), compute magnetic field lines by solving
        $$ \frac{ds}{dl} = \frac{B^s}{B} $$
        $$ \frac{d\theta}{dl} = \frac{B^\theta}{B} $$ 
        $$ \frac{d\zeta}{dl} = \frac{B^\zeta}{B} $$
    Args:
        bField: the toroidal magnetic field. 
        s0: list of s components of initial points. 
        theta0: list of theta components of initial points. 
        zeta0: list of zeta components of initial points. 
        niter: Number of toroidal periods. 
        nstep: Number of intermediate step for one period
    """

    if isinstance(s0, float):
        s0, theta0, zeta0 = np.array([s0]), np.array([theta0]), np.array([zeta0])
    elif isinstance(s0, list):
        s0, theta0, zeta0 = np.array(s0), np.array(theta0), np.array(zeta0)
    assert s0.shape == theta0.shape == zeta0.shape
    if kwargs.get("method") is None:
        kwargs.update({"method": "LSODA"}) 
    if kwargs.get("rtol") is None:
        kwargs.update({"rtol": 1e-10}) 
    print("Change the resolution of the field... ")
    bField.changeResolution(sResolution=sResolution, thetaResolution=thetaResolution, zetaResolution=zetaResolution)
    print("Get the Jacobian and metric of the field... ")
    baseJacobian = bField.getJacobian()
    baseMetric = bField.getMetric()

    from pyoculus.problems import SPECBfield
    pyoculusField = SPECBfield(bField.specData, bField.lvol+1)
    def getB(dLength, point):
        field = pyoculusField.B([point[0], point[1], point[2]]) / bField.interpValue(baseJacobian, point[0], point[1], point[2])
        metric = bField.interpValue(baseMetric, point[0], point[1], point[2])
        # bSupS = field[0]
        # bSupTheta = field[1]
        # bSupZeta = field[2]
        bPow = 0
        for i in range(3):
            for j in range(3): 
                bPow += (field[i]*field[j]*metric[0,i,j]) 
        b = np.power(bPow, 0.5)
        return [field[0]/b, field[1]/b, field[2]/b]

    print("Begin field line tracing: ")
    lines = list()
    for lineIndex in range(len(s0)):           # loop over each field line
        point = [s0[lineIndex], theta0[lineIndex], zeta0[lineIndex]]
        initLength = 0
        deltaLength = oneLength / nstep
        sArr = [s0[lineIndex]]
        thetaArr = [theta0[lineIndex]]
        zetaArr = [zeta0[lineIndex]]
        for j in range(niter):              # loop over each toroidal iteration
            print_progress(lineIndex*niter+j+1, len(s0)*niter)
            for k in range(nstep):          # loop inside one iteration
                sol = solve_ivp(getB, (initLength, initLength+deltaLength), point, **kwargs)            # solve ODEs
                sArr.append(sol.y[0,-1])
                thetaArr.append(sol.y[1,-1])
                zetaArr.append(sol.y[2,-1])
                point = [sArr[-1], thetaArr[-1], zetaArr[-1]]
                initLength += deltaLength
        lines.append(FieldLine.getLine_tracing(bField, nstep, np.array(sArr), np.array(thetaArr), np.array(zetaArr)))
    return lines


if __name__ == "__main__":
    pass
