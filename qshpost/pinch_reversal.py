from mpy import SPECOut
from mpy.specMagneticField import SPECField
import numpy as np
from scipy.integrate import quad
from typing import Tuple


def getPinchReversalPara(spec_file: str, majorRadius: float, minorRadius: float, phiEdge: float) -> Tuple:
    """
    returns:
        pinch parameter, reversal parameter
    """
    speclib = SPECOut(spec_file)
    outerField = SPECField(specData=speclib, lvol=1, sResolution=128, thetaResolution=128, zetaResolution=128)     
    baseJacobian = outerField.getJacobian()
    deltaS = 1e-8
    from pyoculus.problems import SPECBfield
    pyoculusField = SPECBfield(speclib, 2)
    def getBoundaryB(theta, zeta):
        field = pyoculusField.B([1-deltaS, theta, zeta]) / outerField.interpValue(baseJacobian, 1-deltaS, theta, zeta)
        return field[0], field[1], field[2]
    def getB_theta(theta):
        bSupS, bSupTheta, bSupZeta = getBoundaryB(theta, 0)
        return bSupTheta
    def getB_zeta(zeta):
        bSupS, bSupTheta, bSupZeta = getBoundaryB(0, zeta)
        return bSupZeta
    pinchParameter = (
        (quad(getB_theta, 0, 2*np.pi)[0] / (2*np.pi))
        / (phiEdge / (np.pi*minorRadius*minorRadius))
    )
    reversalParameter = (
        (quad(getB_zeta, 0, 2*np.pi)[0] / (2*np.pi))
        / (phiEdge / (np.pi*minorRadius*minorRadius))
    )
    return pinchParameter, reversalParameter


if __name__ == "__main__":
    pass
