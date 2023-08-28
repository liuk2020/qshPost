import numpy as np
from scipy.integrate import solve_ivp 
from mpy.specMagneticField import SPECField, FieldLine
from mpy.specMagneticField import readJacobian, readB
from mpy.misc import print_progress
from typing import List


def traceLine(
    bField: SPECField, 
    s0: np.ndarray, theta0: np.ndarray, zeta0: np.ndarray, 
    niter: int=128, nstep: int=32,
    bMethod: str="calculate", 
    bData: str=None, jacobianData: str=None, 
    printControl: bool=True, **kwargs
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
        bMethod: should be `"calculate"` or "`interpolate`" , the method to get the magnetic field. 
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
    
    if bMethod == "calculate":
        from pyoculus.problems import SPECBfield
        pyoculusField = SPECBfield(bField.specData, bField.lvol+1)
        if jacobianData is None:
            base_Jacobian = bField.getB()
            base_sArr = bField.sArr
            base_thetaArr = bField.thetaArr
            base_zetaArr = bField.zetaArr
        else:
            base_sArr, base_thetaArr, base_zetaArr, base_Jacobian = readJacobian(jacobianData)
    elif bMethod == "interpolate":
        if bData is None:
            base_bSupS, base_bSupTheta, base_bSupZeta = bField.getB()
            base_sArr = bField.sArr
            base_thetaArr = bField.thetaArr
            base_zetaArr = bField.zetaArr
        else:
            base_sArr, base_thetaArr, base_zetaArr, base_bSupS, base_bSupTheta, base_bSupZeta = readB(bData)
    else:
        raise ValueError(
            "`bMethod` should be `calculate` or `interpolate`. "
        )

    def getB_calculate(zeta, s_theta):
        # field = pyoculusField.B_many(s_theta[0], s_theta[1], zeta) / bField.interpValue(base_Jacobian, s_theta[0], s_theta[1], zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        # bSupS = field[0, 0]
        # bSupTheta = field[0, 1]
        # bSupZeta = field[0, 2]
        field = pyoculusField.B([s_theta[0], s_theta[1], zeta]) / bField.interpValue(base_Jacobian, s_theta[0], s_theta[1], zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        bSupS = field[0]
        bSupTheta = field[1]
        bSupZeta = field[2]
        return [bSupS/bSupZeta, bSupTheta/bSupZeta]
    
    def getB_interpolate(zeta, s_theta):
        bSupS = bField.interpValue(baseData=base_bSupS, sValue=s_theta[0], thetaValue=s_theta[1], zetaValue=zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        bSupTheta = bField.interpValue(baseData=base_bSupTheta, sValue=s_theta[0], thetaValue=s_theta[1], zetaValue=zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        bSupZeta = bField.interpValue(baseData=base_bSupZeta, sValue=s_theta[0], thetaValue=s_theta[1], zetaValue=zeta, sArr=base_sArr, thetaArr=base_thetaArr, zetaArr=base_zetaArr)
        return [bSupS/bSupZeta, bSupTheta/bSupZeta]
    
    lines = list()
    nLine = len(s0)
    if printControl:
        print("Begin field-line tracing: ")
    for i in range(nLine):              # loop over each field-line 
        s_theta = [s0[i], theta0[i]]
        zetaStart = zeta0[i]
        dZeta = 2 * np.pi / bField.nfp / nstep
        sArr = [s0[i]]
        thetaArr = [theta0[i]]
        zetaArr = [zeta0[i]]
        for j in range(niter):          # loop over each toroidal iteration
            if printControl:
                print_progress(i*niter+j+1, nLine*niter)
            for k in range(nstep):      # loop inside one iteration
                if bMethod == "calculate":
                    sol = solve_ivp(
                        getB_calculate, 
                        (zetaStart, zetaStart+dZeta), 
                        s_theta, **kwargs
                    )
                elif bMethod == "interpolate":
                    sol = solve_ivp(
                        getB_interpolate, 
                        (zetaStart, zetaStart+dZeta), 
                        s_theta, **kwargs
                    )
                else:
                    raise ValueError("bMethod should be calculate or interpolate. ")
                sArr.append(sol.y[0,-1])
                thetaArr.append(sol.y[1,-1])
                zetaArr.append(zetaStart+dZeta)
                s_theta = [sArr[-1], thetaArr[-1]]
                zetaStart = zetaArr[-1]
        lines.append(FieldLine.getLine_tracing(bField, nstep, np.array(sArr), np.array(thetaArr), np.array(zetaArr)))
    return lines


if __name__ == "__main__":
    pass
