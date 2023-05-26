import h5py
import mpy
import mpy.specMagneticField as specMagneticField
import numpy as np
from scipy.optimize import least_squares
from typing import Tuple


class Axis:

    def __init__(self, initPoint: np.ndarray, xn: np.ndarray, rac: np.ndarray, zas: np.ndarray) -> None:
        """
        R = \sum rac cos(-nv)
        Z = \sum zas sin(-nv)
        """
        self.initPoint = initPoint
        self.xn = xn
        self.rac = rac
        self.zas = zas

    @classmethod
    def readSPECOut(cls, ntor: int, specData: mpy.SPECOut):
        if not specData.input.physics.Istellsym:
            raise ValueError(
                "There is no codes without stellarator symmetry! "
            )
        deltaS = 1e-8
        initPoint = np.array([-1+deltaS, 0, 0])
        xn = specData.output.in_[0: 2*ntor+1]
        rac = specData.output.Rbc[0, 0:2*ntor+1]
        zas = specData.output.Zbs[0, 0:2*ntor+1]
        return cls(initPoint=initPoint, xn=xn, rac=rac, zas=zas)

    @classmethod
    def traceLine(cls, ntor: int, line: specMagneticField.FieldLine, specData: mpy.SPECOut, **kwargs):
        """
        Use the least squares method to fit the axis curve! 
        """
        xn = specData.output.in_[0: 2*ntor+1]
        nums = len(xn)
        if kwargs.get("verbose") is None:
            kwargs.update({"verbose": 1}) 
        def getValue(coeffArr: np.ndarray) -> Tuple[np.ndarray]:
            rac = coeffArr[0: nums]
            zas = coeffArr[nums: 2*nums]
            angleMat = - np.dot(xn.reshape(-1,1), line.zetaArr.reshape(1,-1))
            rArr = (
                np.dot(rac.reshape(1,-1), np.cos(angleMat))
            ).flatten()
            zArr = (
                np.dot(zas.reshape(1,-1), np.sin(angleMat))
            ).flatten()
            return rArr, zArr
        def getErr(coeffArr: np.ndarray) -> np.ndarray:
            rArr, zArr = getValue(coeffArr)
            return np.append(
                np.power(rArr-line.rArr,2), np.power(zArr-line.zArr,2)
            )
        optimizeRes = least_squares(getErr, np.zeros(2*nums), bounds=(-np.inf, np.inf), **kwargs)
        assert optimizeRes.success
        return cls(initPoint=np.array([line.sArr[0],line.thetaArr[0],line.zetaArr[0]]), xn=xn, rac=optimizeRes.x[0:nums], zas=optimizeRes.x[nums:2*nums])

    def getRZ(self, zetaArr: np.ndarray) -> Tuple[np.ndarray]:
        angleMat = - np.dot(self.xn.reshape(-1,1), zetaArr.reshape(1,-1))
        rArr = (
            np.dot(self.rac.reshape(1,-1), np.cos(angleMat))
        ).flatten()
        zArr = (
            np.dot(self.zas.reshape(1,-1), np.sin(angleMat))
        ).flatten()
        return rArr, zArr

    @classmethod
    def readH5(cls, fileName: str):
        with h5py.File(fileName, 'r') as f:
            initPoint = f["initPoint"][:]
            xn = f["xn"][:]
            rac = f["rac"][:]
            zas = f["zas"][:]
        return cls(initPoint=initPoint, xn=xn, rac=rac, zas=zas)

    def writeH5(self, fileName: str):
        with h5py.File(fileName, 'w') as f:
            f.create_dataset("initPoint", data=self.initPoint)
            f.create_dataset("xn", data=self.xn)
            f.create_dataset("rac", data=self.rac)
            f.create_dataset("zas", data=self.zas)


if __name__ == "__main__":
    pass
