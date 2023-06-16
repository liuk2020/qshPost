import numpy as np
from mpy.specMagneticField import FieldLine
from .axis import Axis


def getFirstIota(line: FieldLine) -> float:
    iota = (line.thetaArr[-1]-line.thetaArr[0]) / (line.zetaArr[-1]-line.zetaArr[0])
    return float(iota)


def getSecondIota(line: FieldLine, axis: Axis) -> float:
    
    rAxis, zAxis = axis.getRZ(line.zetaArr)
    rArr = line.rArr - rAxis
    zArr = line.zArr - zAxis
    positionVec = np.concatenate((rArr.reshape(-1,1),zArr.reshape(-1,1)), axis=1)
    nums = len(positionVec)
    normArr, cosAngle = list(), list()
    for i in range(nums):
        normArr.append(np.linalg.norm(positionVec[i]))
    for i in range(nums-1):
        cosAngle.append(np.dot(positionVec[i],positionVec[i+1])/normArr[i]/normArr[i+1])
    thetaArr = np.cumsum(np.arccos(cosAngle))
    iota = thetaArr[-1] / (line.zetaArr[-1]-line.zetaArr[0])
    return float(iota)
