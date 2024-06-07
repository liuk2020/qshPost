import matplotlib.pyplot as plt
from mpy import SPECOut


def plotCase(case: SPECOut, ax=None, toroidalIdx: int=0, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    rr = case.poincare.R[:, :, toroidalIdx]
    zz = case.poincare.Z[:, :, toroidalIdx]
    if kwargs.get("marker") == None:
        kwargs.update({"marker": "."})
    if kwargs.get("c") == None:
        pass
    if kwargs.get("s") == None:
        kwargs.update({"s": 0.3})
    nptrj = rr.shape[0]
    for ii in range(nptrj-1):
        dots = ax.scatter(rr[ii, :], zz[ii, :], **kwargs)
    plt.xlabel(r"$R(\mathrm{m})$", fontsize=20)
    plt.ylabel(r"$Z(\mathrm{m})$", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axis("equal")


def plotQ(case: SPECOut, ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    plt.sca(ax)
    if kwargs.get("marker") == None:
        kwargs.update({"marker": "."})
    xdata = case.poincare.R[:, 0, 0]
    ydata = case.transform.fiota[1, case.poincare.success==1]
    xlabel = r"R"
    ylabel = r"$q$"
    dots = ax.scatter(xdata[1:-1], 1/ydata[1:-1], **kwargs)
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)


if __name__ == "__main__":
    pass
