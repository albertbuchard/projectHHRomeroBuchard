# General utilities


# Project utilities
from utilities import *

# Neuron Simulation Step Current
Istart =10.
Iend =1000.
end=Iend + 50.
Iamp = 8.
HHNeuronSim = HH_Step(I_tstart=Istart, I_tend=Iend, I_amp=Iamp, tend=end, do_plot=True, model='Adapt')

# Get the voltage curve
voltageValues = HHNeuronSim.vm[0]/b2.mV

# x values
sizeX = len(voltageValues)
xvals = list(HHNeuronSim.t)

# convert voltage values to a MathFunction object - using spline smoothing
voltageFunction = MathFunction(voltageValues, "spline", xValues=HHNeuronSim.t)

# find the extremum with some constraints on the voltage and convexity of the second derivative
turningPoints = voltageFunction.FindExtremumInRange(errorValue=voltageFunction.mean / 10000, espaceDistance=0.0005,
                                                    convexity=-1,
                                                    yThreshold=voltageFunction.mean + 40)

print(turningPoints)
plt.figure()
plt.plot(xvals, voltageValues)
for xTurningPoint in turningPoints:
    plt.plot((xTurningPoint, xTurningPoint),
             (voltageFunction(xTurningPoint) - (voltageFunction.sd * 4),
              voltageFunction(xTurningPoint) + (voltageFunction.sd * 4)), 'k-')
plt.show()

# Frequency with respect to time
windowedFrequency = slideCompute(voltageFunction.rangeX, turningPoints, windowSize=0.02, step=0.02)
plt.figure()
plt.plot(windowedFrequency[0], windowedFrequency[2])
plt.show()

HHNeuronSim = None

# FI Curve
plotFICurve(np.arange(0, 40., 0.02), Iend=1000)

# #Computation of Turning points
# voltageValues = HHNeuronSim.vm[0]/b2.mV
# sizeX = len(voltageValues)
# xvals = np.linspace(0,sizeX,sizeX)
#
# # convert voltage values to a MathFunction object - using spline smoothing
# voltageFunction = pu.MathFunction(voltageValues,"spline",)
#
# # find the extremum with some constraints on the voltage and convexity of the second derivative
# turningPoints = voltageFunction.FindExtremumInRange(errorValue=voltageFunction.mean/1000, espaceDistance=5, convexity=-1, yThreshold=voltageFunction.mean+40)
# [windowPositions,numberOfEventsForWindows,frequencyForWindows ]=pu.slideCompute (xRange=[0 ,sizeX], events=turningPoints, step = 100, window = 100, fixedPointOfWindow = "middle")