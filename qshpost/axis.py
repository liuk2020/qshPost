import h5py
import mpy
import mpy.specMagneticField as specMagneticField
import numpy as np
from scipy.integrate import dblquad
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

    def getXYZ(self, zetaArr: np.ndarray) -> Tuple[np.ndarray]:
        rArr, zArr = self.getRZ(zetaArr)
        xArr = np.cos(zetaArr) * rArr
        yArr = np.sin(zetaArr) * rArr
        return xArr, yArr, zArr

    def getWrithe(self) -> float:
        def getRZ(zeta: float) -> Tuple[np.float64]:
            angle = zeta * self.xn.flatten()
            r, z = np.dot(self.rac.flatten(), np.cos(angle)), np.dot(self.zas.flatten(), np.sin(angle))
            return r, z
        def getdRZ(zeta: float) -> Tuple[np.float64]:
            angle = zeta * self.xn.flatten()
            dR = np.dot(self.xn*self.rac.flatten(), np.sin(-angle))
            dZ = np.dot(self.xn*self.zas.flatten(), np.cos(angle))
            return dR, dZ
        def getValue(zeta2, zeta1):
            if zeta1 == zeta2:
                return 0
            r1, z1 = getRZ(zeta1)
            x1, y1 = r1*np.cos(zeta1), r1*np.sin(zeta1)
            r2, z2 = getRZ(zeta2)
            x2, y2 = r2*np.cos(zeta2), r2*np.sin(zeta2)
            dR1, dZ1 = getdRZ(zeta1)
            dX1, dY1 = dR1*np.cos(zeta1)-r1*np.sin(zeta1), dR1*np.sin(zeta1)+r1*np.cos(zeta1)
            dR2, dZ2 = getdRZ(zeta2)
            dX2, dY2 = dR2*np.cos(zeta2)-r2*np.sin(zeta2), dR2*np.sin(zeta2)+r2*np.cos(zeta2)
            dXYZ1 = np.array([dX1, dY1, dZ1])
            dXYZ2 = np.array([dX2, dY2, dZ2])
            deltaXYZ = np.array([x2-x1, y2-y1, z2-z1])
            norm = np.linalg.norm(deltaXYZ)
            return np.linalg.det(np.array([[dX1, dY1, dZ1], [dX2, dY2, dZ2], [x2-x1, y2-y1, z2-z1]])) / np.power(norm,3)        
        value, err = dblquad(getValue,0,2*np.pi,0,2*np.pi)
        print("An estimate of the error is " + str(err) + " ...")
        ans = value/np.pi/4
        print("The writhe of the axis curve is " + str(ans) + " ...")
        return ans

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
