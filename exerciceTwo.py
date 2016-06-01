# Project utilities
from utilities import *

# Neuron Simulation Step Current
Istart =100.
Iend =1000.
end=Iend + 50.
Iamp = 8.

#
# # First model for adaptation only M type channels
# HHNeuronSim = HH_Step(I_tstart=Istart, I_tend=Iend, I_amp=Iamp, tend=end, do_plot=True, model='Monly')
#
# # Get the voltage curve
# voltageValues = HHNeuronSim.vm[0]/b2.mV
#
# # x values
# sizeX = len(voltageValues)
# xvals = list(HHNeuronSim.t)
#
# # convert voltage values to a MathFunction object - using spline smoothing
# voltageFunction = MathFunction(voltageValues, "spline", xValues=HHNeuronSim.t)
#
# # find the extremum with some constraints on the voltage and convexity of the second derivative
# turningPoints = voltageFunction.FindExtremumInRange(errorValue=voltageFunction.mean / 10000, espaceDistance=0.0005,
#                                                     convexity=-1,
#                                                     yThreshold=voltageFunction.mean + 40)
#
# print(turningPoints)
# plt.figure()
# plt.plot(xvals, voltageValues)
# for xTurningPoint in turningPoints:
#     plt.plot((xTurningPoint, xTurningPoint),
#              (voltageFunction(xTurningPoint) - (voltageFunction.sd * 4),
#               voltageFunction(xTurningPoint) + (voltageFunction.sd * 4)), 'k-')
# plt.draw()
#
# # Frequency with respect to time
# windowedFrequency = slideCompute(voltageFunction.rangeX, turningPoints, windowSize=0.05, step=0.05)
# plt.figure()
# plt.plot(windowedFrequency[0],windowedFrequency[2])
# plt.title("Evolution of frequency in time for HH model without M type channel")
# plt.xlabel('t (ms)')
# plt.ylabel('f (Hz)')
# plt.subplot(111).spines['right'].set_color((.8, .8, .8))
# plt.subplot(111).spines['top'].set_color((.8, .8, .8))
# plt.draw()
#
# HHNeuronSim = None
#
# # FI Curve
# plotFICurve(np.arange(0, 40., 5), Iend=1000, model="Monly")
#
# plt.show()
#
#
# # Adaptation with M type channels and AHP
# print("M and AHP analysis")
# HHNeuronSim = HH_Step(I_tstart=Istart, I_tend=Iend, I_amp=Iamp, tend=end, do_plot=True, model='MandAHP')
#
# # Get the voltage curve
# voltageValues = HHNeuronSim.vm[0]/b2.mV
#
# # x values
# sizeX = len(voltageValues)
# xvals = list(HHNeuronSim.t)
#
# # convert voltage values to a MathFunction object - using spline smoothing
# voltageFunction = MathFunction(voltageValues, "spline", xValues=HHNeuronSim.t)
#
# # find the extremum with some constraints on the voltage and convexity of the second derivative
# turningPoints = voltageFunction.FindExtremumInRange(errorValue=voltageFunction.mean / 10000, espaceDistance=0.0005,
#                                                     convexity=-1,
#                                                     yThreshold=voltageFunction.mean + 40)
#
# print(turningPoints)
# plt.figure()
# plt.plot(xvals, voltageValues)
# for xTurningPoint in turningPoints:
#     plt.plot((xTurningPoint, xTurningPoint),
#              (voltageFunction(xTurningPoint) - (voltageFunction.sd * 4),
#               voltageFunction(xTurningPoint) + (voltageFunction.sd * 4)), 'k-')
# plt.draw()
#
# # Frequency with respect to time
# windowedFrequency = slideCompute(voltageFunction.rangeX, turningPoints, windowSize=0.05, step=0.05)
# plt.figure()
# plt.plot(windowedFrequency[0],windowedFrequency[2])
# plt.title("Evolution of frequency in time for HH model with M and AHP type channels")
# plt.xlabel('t (ms)')
# plt.ylabel('f (Hz)')
# plt.subplot(111).spines['right'].set_color((.8, .8, .8))
# plt.subplot(111).spines['top'].set_color((.8, .8, .8))
# plt.draw()
#
# HHNeuronSim = None
#
# # FI Curve
# plotFICurve(np.arange(0, 40., 5), Iend=1000, model="MandAHP")
#
# plt.show()
#

# Adaptation with M type channels, AHP, and Ca channels
print("M, AHP and Ca analysis")
HHNeuronSim = HH_Step(I_tstart=Istart, I_tend=Iend, I_amp=Iamp, tend=end, do_plot=True, model='MAHPandCa')

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
plt.draw()

# Frequency with respect to time
windowedFrequency = slideCompute(voltageFunction.rangeX, turningPoints, windowSize=0.05, step=0.05)
plt.figure()
plt.plot(windowedFrequency[0],windowedFrequency[2])
plt.title("Evolution of frequency in time for HH model with M, AHP and Ca type channels")
plt.xlabel('t (ms)')
plt.ylabel('f (Hz)')
plt.subplot(111).spines['right'].set_color((.8, .8, .8))
plt.subplot(111).spines['top'].set_color((.8, .8, .8))
plt.draw()

HHNeuronSim = None

# FI Curve
plotFICurve(np.arange(0, 40., 5), Iend=1000, model="MAHPandCa")

plt.show()