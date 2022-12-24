import csv
import os
import sys
import math
from matplotlib import animation as anim
from matplotlib import pyplot as plt
from time import time
import scipy.integrate as integrate
import bisect
import numpy as npy

'''
class circleDrawer(object):
    @staticmethod
    def periodicfunc(begin, end, func, x):
        # todo use mode function?
        if x >= begin and x <= end:
            return func(x)
        elif x > end:
            xNew = x - (end - begin)
            return FourierSeries.periodicfunc(begin, end, func, xNew)
        elif x < begin:
            xNew = x + (end - begin)
            return FourierSeries.periodicfunc(begin, end, func, xNew)

    def __init__(self, an, bn, x=0.0, y=0.0):
        self.an_ = an
        self.bn = bn
        self.x = x
        self.y = y

 def animate(self, i, xVec:list, yVec:list):
        r = sqrt(())

        for i in range(0, self.numP):
            t_ =
            xVec[i] = self.x + npy.cos(t_) * r
            yVec[i] = self.y + npy.sin(t_) * r
            anim.'''


class FourierSeries(object):
    @staticmethod
    def isFloat(val):
        try:
            float(val)
            return True
        except ValueError:
            return False

    def __init__(self, fileName,  numPt, integrationRule):
        self._numPt = numPt
        self._shapeFile = fileName
        self._integrationRule = integrationRule

    def LoadShape(self, fileName):
        folder = os.path.dirname(os.path.realpath(fileName))
        fullPath = (os.path.join(folder, os.path.normpath(fileName))).replace('\\','/')
        normalizedFullPath = os.path.normpath(fullPath)

        if os.path.isfile(normalizedFullPath) is False:
            raise FileExistsError("LoadShap: file: {} not exists".format(normalizedFullPath))

        xVec = []
        yVec = []
        distVec = []

        with open(normalizedFullPath, 'r', encoding='utf-8-sig') as data:
            rows = csv.DictReader(data)
            i = 0

            for row in rows:
                if len(row) != 2:
                    raise RuntimeError("invalid input file")

                x = row['x'].strip()
                y = row['y'].strip()

                if FourierSeries.isFloat(x) is False or FourierSeries.isFloat(y) is False:
                    raise ValueError("invalid values in the input file: {} {}".format(x, y))

                xVec.append(float(int(x)))
                yVec.append(-float(int(y)))
                dist = 0 if i == 0 else distVec[i-1] + math.sqrt((xVec[i]-xVec[i-1])**2+(yVec[i]-yVec[i-1])**2)
                distVec.append(float(dist))
                i += 1

        data.close()

        return xVec, yVec, distVec

    @staticmethod
    def Integrate(func, begin, end, args):

        i, harmonic, xVec, yVec, period, num= args

        #num = 250

        val = 0.0

        delta = (end - begin)/(num-1)
        for index in range(1, num):
            val += delta * (func(index * delta, i, harmonic, xVec, yVec, period) +
                    func((index-1) * delta, i, harmonic, xVec, yVec, period)) * 0.5

        return val

    @staticmethod
    def integrand(x, i, harmonicFunc, xVec, yVec, T):
        #todo: use numpy.interp?

        if x < xVec[0] or x > xVec[-1]:
            return 0

        pos = bisect.bisect_left(xVec, x)
        yVal = yVec[pos-1] + (yVec[pos]-yVec[pos-1])/(xVec[pos]-xVec[pos-1])*(x-xVec[pos-1])

        h = harmonicFunc(2.0 * i * npy.pi * x / T)

        return yVal if i == 0 else yVal * h

    def findFTC(self, points, n=0):

        if n == 0:
            n = len(points[0])

        # period
        begin = points[2][0]
        end = points[2][-1]

        period = end - begin

        # check the integrand
        #num = 2000
        #tempVec = [None] * num
        #tempX = [None] * num

        #for i in range(0, num):
        #    tempX[i] = i* period/(num-1.0)
        #    tempVec[i] = FourierSeries.integrand(tempX[i], 0, npy.cos, points[2], points[0], period), self._numPt

        #plt.plot(tempX, tempVec)
        #plt.show()

        invP = 1.0 / period

        # x Coord
        # const term

        a0_x = invP * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                             args=(0, npy.cos, points[2], points[0], period, self._numPt))
        #a0_x = invP * integrate.quad(FourierSeries.integrand, begin, end,
        #                                    args=(0, npy.cos, points[2], points[0], period, self._numPt))

        print("Coefficients for X Component:")
        print("A0: {}".format(a0_x))
        A_x = npy.zeros((n-1))
        B_x = npy.zeros((n-1))


        for i in range(1, n):
            A_x[i - 1] = invP * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                                       args=(i, npy.cos, points[2], points[0], period, self._numPt))
            B_x[i - 1] = invP * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                                       args=(i, npy.sin, points[2], points[0], period, self._numPt))
            print("{}, {}".format(A_x[i-1], B_x[i-1]))

        # y coord

        print("Coefficients for Y Component:")
        a0_y = invP  * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                            args=(0, npy.cos, points[2], points[1], period, self._numPt))

        print("A0: {}".format(a0_y))
        A_y = npy.zeros((n-1))
        B_y = npy.zeros((n-1))

        for i in range(1, n):
            A_y[i - 1] = invP * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                                       args=(i, npy.cos, points[2], points[1], period, self._numPt))
            B_y[i - 1] = invP * FourierSeries.Integrate(FourierSeries.integrand, begin, end,
                                                       args=(i, npy.sin, points[2], points[1], period, self._numPt))

            print("{}, {}".format(A_y[i - 1], B_y[i - 1]))

        return a0_x*0.5, A_x, B_x, a0_y*0.5, A_y, B_y

    def Invert(self, num, period, args, N_max=0):

        a0_x, An_x, Bn_x, a0_y, An_y, Bn_y = args

        if N_max == 0:
            N_max = len(An_x)

        delta = period/(num-1)
        xVec = [None] * num
        yVec = [None] * num


        for i in range(0, num):
            xVal = a0_x
            yVal = a0_y
            t = delta * i
            for n in range(1, N_max):
                t_ = n * 2.0 * npy.pi * t / period
                c_h = npy.cos(t_)
                s_h = npy.sin(t_)

                xVal += An_x[n-1] * c_h + Bn_x[n-1] * s_h
                yVal += An_y[n-1] * c_h + Bn_y[n-1] * s_h

            xVec[i] = xVal
            yVec[i] = yVal


        return xVec, yVec;


    def Invoke(self):
        points = [[], [], []]
        points[0], points[1], points[2] = self.LoadShape(self._shapeFile)

        # draw x,y
        plt.title("Original Shape--IB Logo")
        plt.plot(points[0], points[1])
        plt.show()

        # draw d,x
        plt.title("Parameterized X and Y Curve")
        plt.plot(points[2], points[0])
        # draw d,y
        plt.plot(points[2], points[1])
        plt.show()

        # find Fourier coefficients
        coefficients = (0.0, [], [], 0.0)
        coefficients = self.findFTC(points, 128)

        # inverse
        period = points[2][-1] - points[2][0]
        print("period: {}".format(period))

        xVec, yVec= self.Invert(256, period, coefficients)

        # draw recovered shape

        plt.title("Recovered Shape--IB Logo")
        plt.plot(xVec, yVec)
        plt.show()

        # draw animation of circle plot

        #t0 = time()
        #animate(0)
        #t1 = time()
        #interval = numPt * dt - (t1-t0)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    numArgv = len(sys.argv)

    numPStr = sys.argv[1] if numArgv >= 2 else "250"
    numPt = int(numPStr)

    if numPt < 140 or numPt > 500:
        sys.exit("warning: bad number of points: {}".format(numPt))

    path_to_shape =''

    if numArgv >= 3:
        fileN = sys.argv[2]
        bundle_dir = os.getcwd()
        print("bundle_dir: {}".format(bundle_dir))
        path_to_shape = os.path.abspath(os.path.join(bundle_dir, fileN))
        print("path_to_shape: {}".format(path_to_shape), flush=True)
    else:
        bundle_dir = getattr(sys, '_MEIPASS', os.path.abspath(os.path.dirname(__file__)))
        print("no file specified. bundle_dir: {}".format(bundle_dir))
        path_to_shape = os.path.abspath(os.path.join(bundle_dir, 'newIB.txt'))
        print("path_to_shape: {}".format(path_to_shape), flush=True)

    fs = FourierSeries(path_to_shape, numPt, 'trapezoidal')
    fs.Invoke()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
