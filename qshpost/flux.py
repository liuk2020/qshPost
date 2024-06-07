import numpy as np 
import mpy.specMagneticField as specMagneticField 
# import mpy.fitting as fitting
from scipy.integrate import quad
from pyoculus.problems import SPECBfield
# from mpy.specMagneticField import FieldLine
from .crossSurface import FirstCrossSurface, SecondCrossSurface


def getFirstFlux(bField: specMagneticField.SPECField, crossSurf: FirstCrossSurface) -> float: 
    
    pyoculusField = SPECBfield(bField.specData, bField.lvol+1)
    def getPotential(theta):
        return pyoculusField.vectorPotential([crossSurf.getS(np.array([theta])),np.array([theta]),np.zeros(1)])[0]
    
    res = quad(getPotential, 0, 2*np.pi)
    return res[0]


def getSecondFlux(bField: specMagneticField.SPECField, crossSurf: SecondCrossSurface) -> float: 

    pyoculusField = SPECBfield(bField.specData, bField.lvol+1)
    def getPotential(label):
        return pyoculusField.vectorPotential([crossSurf.getS(np.array([label])),crossSurf.getTheta(np.array([label])),np.zeros(1)])[0]
    
    def dTheta(label):
        deltaLabel = 1e-8
        return (crossSurf.getTheta(np.array([label+deltaLabel]))-crossSurf.getTheta(np.array([label-deltaLabel]))) / deltaLabel / 2

    def getValue(label):
        return dTheta(label) * getPotential(label)
    
    res = quad(getValue, 0, 2*np.pi)
    return res[0]


if __name__ == "__main__": 
    pass
