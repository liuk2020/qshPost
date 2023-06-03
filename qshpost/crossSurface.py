import numpy as np 
import mpy.fitting as fitting
from scipy.optimize import least_squares
from mpy.specMagneticField import FieldLine
from typing import Tuple


class CrossSurface:

    def __init__(self, xm: np.ndarray, rc: np.ndarray, rs: np.ndarray, zc: np.ndarray, zs: np.ndarray) -> None: 
        self.xm = xm
        self.rc = rc
        self.rs = rs
        self.zc = zc
        self.zs = zs

    @classmethod
    def fitFirstCross(cls, line: FieldLine, mpol: int, **kwargs):
        rArr, zArr, thetaArr = list(), list(), list()
        for i in range(len(line.rArr)):
            if i % line.nZeta == 0:
                rArr.append(line.rArr[i])
                zArr.append(line.zArr[i])
                thetaArr.append(line.thetaArr[i])
        xm, rc = fitPeriodicR(np.array(thetaArr), np.array(rArr), mpol, **kwargs) 
        xm, zs = fitPeriodicZ(np.array(thetaArr), np.array(zArr), mpol, **kwargs) 
        return cls(xm=xm, rc=rc, rs=np.zeros(len(xm)), zc=np.zeros(len(xm)), zs=zs)

    @classmethod
    def fitSecondCross(cls, line: FieldLine, mpol: int, **kwargs):
        rArr, zArr = list(), list()
        for i in range(len(line.rArr)):
            if i % line.nZeta == 0 and line.zArr[i] >= 0:
                rArr.append(line.rArr[i])
                zArr.append(line.zArr[i])
        rArr = np.array(rArr)
        zArr = np.array(zArr)
        nums = rArr.size
        pointState = [False for _i in range(nums)]
        pointIndex = 0
        pointState[pointIndex]= True
        _rArr, _zArr = [rArr[pointIndex]], [zArr[pointIndex]]
        for i in range(nums):
            minDistance = 1e10
            _index = pointIndex
            for j in range(nums):
                if pointState[j] == False:
                    distance = np.power(rArr[pointIndex]-rArr[j],2) +  np.power(zArr[pointIndex]-zArr[j],2)
                    if distance < minDistance:
                        minDistance = distance
                        _index = j
            pointIndex = _index
            pointState[pointIndex]= True
            _rArr.append(rArr[pointIndex])
            _zArr.append(zArr[pointIndex])
        for i in range(nums):
            _rArr.append(_rArr[nums-1-i])
            _zArr.append(-_zArr[nums-1-i])
        disArr = np.zeros(len(_rArr))
        for i in range(len(_rArr)-1):
            disArr[i+1] = np.power(np.power(_rArr[i+1]-_rArr[i],2)+np.power(_zArr[i+1]-_zArr[i],2), 0.5)
        sumDis = np.cumsum(disArr)
        _thetaArr = 2*np.pi*sumDis/sumDis[-1]
        xm, rc = fitPeriodicR(_thetaArr, np.array(_rArr), mpol)
        xm, zs = fitPeriodicZ(_thetaArr, np.array(_zArr), mpol) 
        return cls(xm=xm, rc=rc, rs=np.zeros(len(xm)), zc=np.zeros(len(xm)), zs=zs)

    @classmethod 
    def fitFirstCross_old(cls, line: FieldLine, mpol: int):
        rArr, zArr, thetaArr = list(), list(), list()
        for i in range(len(line.rArr)):
            if i % line.nZeta == 0:
                rArr.append(line.rArr[i])
                zArr.append(line.zArr[i])
                thetaArr.append(line.thetaArr[i])
        xm, rs, rc = fitting.fitPeriodicCurve(np.array(thetaArr), np.array(rArr), mpol)
        xm, zs, zc = fitting.fitPeriodicCurve(np.array(thetaArr), np.array(zArr), mpol)
        return cls(xm=xm, rc=rc, rs=rs, zc=zc, zs=zs)

    def getRZ(self, theta: np.ndarray) -> Tuple[np.ndarray]: 
        angleMat = np.dot(self.xm.reshape(-1,1), theta.reshape(1,-1))
        rArr = (
            np.dot(self.rc.reshape(1,-1), np.cos(angleMat)) + 
            np.dot(self.rs.reshape(1,-1), np.sin(angleMat))
        ).flatten()
        zArr = (
            np.dot(self.zc.reshape(1,-1), np.cos(angleMat)) + 
            np.dot(self.zs.reshape(1,-1), np.sin(angleMat))
        ).flatten()
        return rArr, zArr


def fitPeriodicR(thetaArr: np.ndarray, rArr: np.ndarray, mpol: int, debug: bool=False, **kwargs) -> Tuple[np.ndarray]:
    """
    Use the least squares method to fit the periodic curve! 
        r = \sum(rc*cos(xm*theta))
    return:
        xm, rc 
    """

    assert thetaArr.shape == rArr.shape
    assert (mpol+1) < thetaArr.size
    # `verbose = 1`: display a termination report.(0: work silently; 2: display progress during iterations (not supported by 'lm' method))
    if kwargs.get("verbose") is None:
        kwargs.update({"verbose": 1}) 
    xm = np.arange(mpol+1)
    
    def getS(rc: np.ndarray, thetaArr: np.ndarray) -> np.ndarray:
        angleMat = np.dot(xm.reshape(-1,1), thetaArr.reshape(1,-1))
        return (
            np.dot(rc.reshape(1,-1), np.cos(angleMat))
        ).flatten()
    
    def getErr(rc: np.ndarray, angleArr: np.ndarray, s: np.ndarray) -> np.ndarray:
        return s - getS(rc, angleArr)
    
    optimizeRes = least_squares(getErr, np.zeros(mpol+1), args=(thetaArr, rArr), bounds=(-np.inf, np.inf), **kwargs)

    if debug:
        return optimizeRes
    else:
        assert optimizeRes.success
        return xm, optimizeRes.x[:]


def fitPeriodicZ(thetaArr: np.ndarray, zArr: np.ndarray, mpol: int, debug: bool=False, **kwargs) -> Tuple[np.ndarray]:
    """
    Use the least squares method to fit the periodic curve! 
        z = \sum(rc*cos(xm*theta))
    return:
        xm, zs 
    """

    assert thetaArr.shape == zArr.shape
    assert (mpol+1) < thetaArr.size
    # `verbose = 1`: display a termination report.(0: work silently; 2: display progress during iterations (not supported by 'lm' method))
    if kwargs.get("verbose") is None:
        kwargs.update({"verbose": 1}) 
    xm = np.arange(mpol+1)
    
    def getS(zs: np.ndarray, thetaArr: np.ndarray) -> np.ndarray:
        angleMat = np.dot(xm.reshape(-1,1), thetaArr.reshape(1,-1))
        return (
            np.dot(zs.reshape(1,-1), np.sin(angleMat))
        ).flatten()
    
    def getErr(zs: np.ndarray, angleArr: np.ndarray, s: np.ndarray) -> np.ndarray:
        return s - getS(zs, angleArr)
    
    optimizeRes = least_squares(getErr, np.zeros(mpol+1), args=(thetaArr, zArr), bounds=(-np.inf, np.inf), **kwargs)

    if debug:
        return optimizeRes
    else:
        assert optimizeRes.success
        return xm, optimizeRes.x[:]


if __name__ == "__main__":
    pass
