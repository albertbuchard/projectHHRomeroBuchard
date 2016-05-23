import numpy as np
import matplotlib.pyplot as plt
import brian2 as b2
from scipy import interpolate
import random as rand   



# Math functions
def plotFICurve(iAmpList, model = "HH", end = 270.0, Istart= 10.0, Iend = 260.0, Iamp = 8.0):

    # set window length to the duration of the step current, scale it to seconds
    windowSize = (Iend - Istart)/1000

    # setup output frequency list
    frequencyList = []

    # loop over the current parameters
    for iAmp in iAmpList:
        # Do the ramp and get the voltage and time points
        neuron = None
        neuron = HH_Step(I_amp=iAmp, I_tstart=Istart, I_tend=Iend, tend=end, do_plot=False, model=model)

        voltageValues = neuron.vm[0]
        xvals = list(neuron.t)

        # convert voltage values to a MathFunction object - using spline smoothing
        voltageFunction = MathFunction(voltageValues, "spline", xValues=xvals)

        # find the extremum with some constraints on the voltage and convexity of the second derivative
        turningPoints = voltageFunction.FindExtremumInRange(errorValue=voltageFunction.mean / 10000,
                                                            espaceDistance=0.0005, convexity=-1,
                                                            yThreshold=voltageFunction.mean + voltageFunction.sd * 1.5)

        # store average frequency
        frequencyList.append(len(turningPoints) / windowSize)

    plt.figure()
    plt.plot(iAmpList, frequencyList)
    plt.show()

def slideCompute (xRange, events, step = None, windowSize = 99, fixedPointOfWindow = "middle") :

    # Function returning a 3 dimensional list
    # [[0]] windowPosition : default to the middle of the window
    # [[1]] The number of event falling into each windows
    # [[2]] The frequency for the window

    # setup output variables
    windowPositions = []
    numberOfEventsForWindows = []
    frequencyForWindows = []

    # define a step if none is given
    if step is None:
        steps = (xRange[1]-xRange[0])/1000

    # start loop over the windows
    for xBegin in np.arange(xRange[0], xRange[1], step) :
        # casting to float
        xBegin = float(xBegin)

        # set the end position of the window
        xEnd = xBegin + windowSize
        if xEnd>xRange[1] :
            xEnd = xRange[1]

        # find the fixed point defining the window for the final returned list
        if fixedPointOfWindow == "middle" :
            xPositionOfWindow = (xEnd + xBegin) / 2
        elif fixedPointOfWindow == "beginning":
            xPositionOfWindow = xBegin

        # check for the number of events in the window
        numberOfEvents = sum(xBegin < x < xEnd for x in events)

        # append output values for window
        windowPositions.append(xPositionOfWindow)
        numberOfEventsForWindows.append(numberOfEvents)
        frequencyForWindows.append(float(numberOfEvents)/float(xEnd-xBegin))

    return [windowPositions, numberOfEventsForWindows, frequencyForWindows]

class MathFunction:
    yValues = []
    xValues = []
    rangeX = []
    mean = None
    median = None
    sd = None
    smoothing = ''
    __smoothedFunction = None
    __splineObject = None
    __validSmoothing = ["none", "spline", "linear", "cubic"]

    def __init__(self, yValues, smoothing, xValues=None):
        yValues = [float(i) for i in yValues]
        self.yValues = yValues
        self.mean = np.mean(yValues)
        self.median = np.median(yValues)
        self.sd = np.std(yValues)

        if xValues is None:
            self.xValues = range(0, len(yValues))
        else:
            xValues = list(xValues)
            xValues = [float(i) for i in xValues]
            if len(xValues) != len(yValues):
                lenList = [len(xValues), len(yValues)]
                self.xValues = xValues[0:min(lenList)]
            else:
                self.xValues = xValues

        self.rangeX = [np.nanmin(self.xValues), np.nanmax(self.xValues)]


        if smoothing in self.__validSmoothing:
            self.smoothing = smoothing

            if smoothing == "spline":
                self.__splineObject = interpolate.splrep(self.xValues, self.yValues)
            else:
                self.__smoothedFunction = interp1d(self.xValues, self.yValues, kind=smoothing)
        else:
            raise Exception(smoothing + " is not a valid smoothing algorithm")

    def __call__(self, x, der=0):
       return self.ForX(x,der)

    def ForX(self, x, der=0):

        # Class method returning the value of the function (der=0) or its derivative (der>0) (for spline smoothing)

        if der!=0:
            if self.__splineObject == None:
                raise Exception("Derivative only implemented for spline smoothing")

        if x<self.xValues[0]:
            return self(self.xValues[0],der)
        # Depending on the interpolation we return different values
        if (self.smoothing=="none"):
            if x in self.xValues:
                return self.yValue[self.xValues.index(x)]
            else:
                previousValue = self.xValues[0]
                closestXIndex = 0
                for i in self.xValues:
                    if sign(self.xValues[i]-x) != sign(previousValue-x):
                        closestXIndex = i
                return self.xValues[closestXIndex]
        elif (self.__smoothedFunction!=None):
            return self.__smoothedFunction(x)
        elif (self.__splineObject != None):
            tck = self.__splineObject
            return(interpolate.splev(x, tck, der=der))
        else:
            raise Exception(self.smoothing + " is an unknown type of smoothing")

    def FindExtremumInRange(self, xRange = None, errorValue = None, espaceDistance = 0, yThreshold = None, convexity=None, step= None):
        if step is None:
            step = self.xValues[1]-self.xValues[0]

        # take the range of all values if not otherwise specified
        if xRange == None:
            xRange = np.arange(self.xValues[0], self.xValues[len(self.xValues)-1],step)

        # Set the errorValue for turning point detection at mean/10000 if not otherwise specified
        if errorValue == None:
            errorValue = self.mean/10000

        # initialize algorithm values
        previousDy = self(self.xValues[0],1)
        turningPoints = []
        currentTurningPoint = None

        previousX = xRange[0]
        for x in xRange[1:]:
            # previous values
            previousY = self(previousX, 0)
            previousDY = self(previousX, 1)
            previousDdY = self(previousX, 2)

            # current values
            y = self(x,0)
            dy = self(x, 1)
            ddy = self(x,2)


            if currentTurningPoint != None:
                # If a turning point is detected wait until a certain escapeDistance
                if (x - currentTurningPoint) > espaceDistance:
                    print("Added turning point at " + str(currentTurningPoint))
                    turningPoints.append(currentTurningPoint)
                    currentTurningPoint = None

            #check if first derivation passes through 0 defining a turning point of the function
            if (np.sign(dy)!=np.sign(previousDY))or(np.abs(dy)<=errorValue):
                # check if a turning point is currently being checked
                if currentTurningPoint != None:
                    # Another turning point was found before escape distance
                    # so current turning point is discarded (protect against noise in an overall constant function)
                    if (x-currentTurningPoint)<=espaceDistance:
                        print("Current turning point at " + str(currentTurningPoint) + " did not escape")
                        # check if the turning point respect the constraints
                        # if y(x) < yThreshold do not count this turning point as valid and continue the loop
                        if ((yThreshold != None) & (y < yThreshold)):
                            currentTurningPoint = None
                            continue

                        # if the convexity (second derivative) of the turning point is not as specified continue the loop
                        if ((convexity != None) & (np.sign(ddy) != convexity)):
                            currentTurningPoint = None
                            continue

                        currentTurningPoint = x

                # if  not add the turning point for testing
                else:
                    # check if the turning point respect the constraints
                    # if y(x) < yThreshold do not count this turning point as valid and continue the loop
                    if ((yThreshold!=None)&(y<yThreshold)):
                        continue

                    # if the convexity (second derivative) of the turning point is not as specified continue the loop
                    if ((convexity!=None)&(np.sign(ddy)!=convexity)):
                        continue

                    currentTurningPoint=x

            previousX = x


        return turningPoints


# Plot functions
def plot_data_Adapt(rec, title=None):
    """Plots a TimedArray for values I and v

    Args:
        rec (TimedArray): the data to plot
        title (string, optional): plot title to display
    """
    plt.figure()
    plt.subplot(311)
    plt.plot(rec.t / b2.ms, rec.vm[0] / b2.mV, lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')
    plt.grid()
    traceall = np.append(rec.m[0], [rec.n[0], rec.h[0], rec.q[0], rec.w[0], rec.sinf[0]])

    plt.subplot(312)

    plt.plot(rec.t / b2.ms, rec.m[0] / np.max(rec.m[0]), 'black', lw=2)
    plt.plot(rec.t / b2.ms, rec.n[0] / np.max(rec.n[0]), 'blue', lw=2)
    plt.plot(rec.t / b2.ms, rec.h[0] / np.max(rec.h[0]), 'red', lw=2)
    plt.plot(rec.t / b2.ms, rec.q[0] / np.max(rec.q[0]), 'green', lw=2)
    plt.plot(rec.t / b2.ms, rec.w[0] / np.max(rec.w[0]), 'magenta', lw=2)
    plt.plot(rec.t / b2.ms, rec.sinf[0] / np.max(rec.sinf[0]), 'cyan', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('act./inact.')
    plt.legend(('m', 'n', 'h', 'q', 'w', 'sinf'))
    plt.grid()
    plt.axis((
        0,
        np.max(rec.t / b2.ms),
        0,
        1.1))

    #    plt.subplot(413)
    #    plt.plot(rec.t/b2.ms, rec.I_e[0]/b2.uamp, lw=2)
    #    plt.axis((
    #        0,
    #        np.max(rec.t/b2.ms),
    #        min(rec.I_e[0]/b2.uamp)*1.1,
    #        max(rec.I_e[0]/b2.uamp)*1.1
    #    ))
    #    plt.xlabel('t (ms)')
    #    plt.ylabel('I (uA)')
    #    plt.grid()
    #
    plt.subplot(313)
    plt.plot(rec.t / b2.ms, -rec.I_K[0] / b2.mamp, 'blue', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_Na[0] / b2.mamp, 'red', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_Ca[0] / b2.mamp, 'cyan', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_M[0] / b2.mamp, 'magenta', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_AHP[0] / b2.mamp, 'green', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('I  (mA)')
    plt.legend(('I K+', 'I Na+', 'I Ca+2', 'I M-type', 'I AHP'))
    plt.grid()

    plt.figure()

    plt.subplot(411)
    plt.plot(rec.t / b2.ms, rec.vm[0] / b2.mV, lw=2)
    plt.xlabel('t [ms]')
    plt.ylabel('v [mV]')
    plt.grid()

    plt.subplot(412)
    plt.plot(rec.t / b2.ms, -rec.I_K[0] / b2.mamp, 'green', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_Na[0] / b2.mamp, 'blue', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_Ca[0] / b2.mamp, 'cyan', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('I  (mA)')
    plt.legend(('I K+', 'I Na+', 'I Ca+2'))
    plt.grid()

    plt.subplot(413)

    plt.plot(rec.t / b2.ms, -rec.I_M[0] / b2.uamp, 'magenta', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_AHP[0] / b2.uamp, 'green', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('I  (uA)')
    plt.legend(('I M-type', 'I AHP'))
    plt.grid()

    plt.subplot(414)

    plt.plot(rec.t / b2.ms, rec.Ca[0], '0.8', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('Ca')
    plt.legend('Ca')
    plt.grid()

    if title is not None:
        plt.suptitle(title)

    plt.show()

def plot_data_S(rec, title=None):
    """Plots a TimedArray for values I and v
    Args:
        rec (TimedArray): the data to plot
        title (string, optional): plot title to display
    """
    plt.figure()
    plt.subplot(411)
    plt.plot(rec.t / b2.ms, rec.vm[0] / b2.mV, lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('V (mV)')
    plt.grid()
    traceall = np.append(rec.m[0], [rec.n[0], rec.h[0]])
    nrmfactor = np.max(traceall) / b2.mV

    plt.subplot(412)

    plt.plot(rec.t / b2.ms, rec.m[0] / nrmfactor / b2.mV, 'black', lw=2)
    plt.plot(rec.t / b2.ms, rec.n[0] / nrmfactor / b2.mV, 'blue', lw=2)
    plt.plot(rec.t / b2.ms, rec.h[0] / nrmfactor / b2.mV, 'red', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('act./inact.')
    plt.legend(('m', 'n', 'h'))
    plt.grid()
    plt.axis((
        0,
        np.max(rec.t / b2.ms),
        0,
        1.1))

    plt.subplot(413)
    plt.plot(rec.t / b2.ms, rec.I_e[0] / b2.uamp, lw=2)
    plt.axis((
        0,
        np.max(rec.t / b2.ms),
        min(rec.I_e[0] / b2.uamp) * 1.1,
        max(rec.I_e[0] / b2.uamp) * 1.1
    ))
    plt.xlabel('t (ms)')
    plt.ylabel('I (uA)')
    plt.grid()

    plt.subplot(414)
    plt.plot(rec.t / b2.ms, -rec.I_K[0] / b2.mamp, 'blue', lw=2)
    plt.plot(rec.t / b2.ms, -rec.I_Na[0] / b2.mamp, 'red', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('I  (mA)')
    plt.legend(('I K+', 'I Na+'))
    plt.grid()

    if title is not None:
        plt.suptitle(title)

    plt.show()

def plot_singleplots(rec):
    plt.figure()
    plt.plot(rec.t / b2.ms, rec.vm[0] / b2.mV, lw=2)

    plt.xlabel('t [ms]')
    plt.ylabel('v [mV]')
    plt.grid()
    plt.title('Membrane Potential [mV]')
    # find max of activation and inactivation variables
    traceall = np.append(rec.m[0], [rec.n[0], rec.h[0]])
    nrmfactor = np.max(traceall) / b2.mV

    plt.figure()

    plt.plot(rec.t / b2.ms, rec.m[0] / nrmfactor / b2.mV, 'black', lw=2)
    plt.plot(rec.t / b2.ms, rec.n[0] / nrmfactor / b2.mV, 'blue', lw=2)
    plt.plot(rec.t / b2.ms, rec.h[0] / nrmfactor / b2.mV, 'red', lw=2)
    plt.xlabel('t (ms)')
    plt.ylabel('act./inact.')
    plt.legend(('m', 'n', 'h'))
    plt.title('Activation Variables')
    plt.grid()


# Neuron definition functions
def HH_Neuron(curr, simtime):

    """Simple Hodgkin-Huxley neuron implemented in Brian2.

    Args:
        curr (TimedArray): Input current injected into the HH neuron
        simtime (float): Simulation time [seconds]

    Returns:
        StateMonitor: Brian2 StateMonitor with recorded fields
        ['vm', 'I_e', 'm', 'n', 'h']
    """

    # neuron parameters
    El = (-67) * b2.mV
    EK = (-100) * b2.mV
    ENa = (50) * b2.mV
    gl = 0.18 * b2.msiemens
    gK = 80 * b2.msiemens
    gNa = 100 * b2.msiemens
    C = 1 * b2.ufarad
    a = 0 * b2.mV

    # forming HH model with differential equations
    eqs = '''

    I_e = curr(t) : amp
    membrane_Im = I_e + gNa*m**3*h*(ENa-vm) + \
        gl*(El-vm) + gK*n**4*(EK-vm)  : amp

    I_Na=gNa*m**3*h*(ENa-vm) : amp
    dm/dt = alpham*(1-m)-betam*m : 1
    alpham = .32*(54*mV+vm)/(1-exp(-(vm+54*mV)/(4*mV)))/mV/ms : Hz
    betam = 0.28*((vm+27*mV)/mV) / (exp((vm+27*mV)/(5*mV))-1)/ms : Hz
    alphah = .128*exp(-(vm+50*mV)/(18*mV))/ms    : Hz
    betah = 4./(1+exp(-(vm+27*mV)/(5*mV)))/ms : Hz
    dh/dt = alphah*(1-h)-betah*h : 1

    I_K= gK*n**4*(EK-vm) : amp
    alphan = .032*(52*mV+vm)/(1-exp(-(vm+52*mV)/(5*mV)))/mV/ms : Hz
    betan = .5*exp(-(vm+57*mV)/(40*mV))/ms : Hz
    dn/dt = alphan*(1-n)-betan*n : 1

    dvm/dt = membrane_Im/C : volt


    taum=1/(alpham+betam)/ms : 1
    taun=1/(alphan+betan)/ms : 1
    tauh=1/(alphah+betah)/ms : 1

    minf=alpham/(alpham+betam) :1
    ninf=alphan/(alphan+betan) :1
    hinf=alphah/(alphah+betah) :1
    '''
    # (Ca/(30+Ca)) : 1
    neuron = b2.NeuronGroup(1, eqs, method='exponential_euler')

    # parameter initialization
    neuron.vm = -67.5 * b2.mV
    neuron.m = 0.0529324852572
    neuron.h = 0.95
    neuron.n = 0.05
    # tracking parameters
    rec = b2.StateMonitor(neuron,
                          ['vm', 'I_e', 'I_K', 'I_Na', 'm', 'n', 'h', 'minf', 'ninf', 'hinf', 'taum', 'tauh',
                           'taun'], record=True)

    # running the simulation
    b2.run(simtime)
    return rec

def HH_Neuron_Adapt(curr, simtime):

    # neuron parameters
    El = (-67) * b2.mV
    EK = (-100) * b2.mV
    ENa = (50) * b2.mV
    ECa = 120 * b2.mV
    EM = -100 * b2.mV
    EAHP = -100 * b2.mV
    gl = 0.2 * b2.msiemens
    gK = 80 * b2.msiemens
    gNa = 100 * b2.msiemens
    gCa = 5 * b2.msiemens
    gM = 3 * b2.msiemens
    gAHP = 3 * b2.msiemens
    C = 1 * b2.ufarad
    a = 0 * b2.mV

    # \+gM*w*(EM-vm)
    # + gAHP*q*(EAHP-vm) +gCa*sinf*(ECa-vm)
    # forming HH model with differential equations with additional Channels Ca, AHP, M type
    eqs = '''

    I_e = curr(t) : amp
    membrane_Im = I_e + gNa*m**3*h*(ENa-vm)+gM*w*(EM-vm)  \
         +gl*(El-vm) + gK*n**4*(EK-vm): amp

    I_Na=gNa*m**3*h*(ENa-vm) : amp
    dm/dt = alpham*(1-m)-betam*m : 1
    alpham = .32*(54*mV+vm)/(1-exp(-(vm+54*mV)/(4*mV)))/mV/ms : Hz
    betam = 0.28*((vm+27*mV)/mV) / (exp((vm+27*mV)/(5*mV))-1)/ms : Hz
    alphah = .128*exp(-(vm+50*mV)/(18*mV))/ms    : Hz
    betah = 4./(1+exp(-(vm+27*mV)/(5*mV)))/ms : Hz
    dh/dt = alphah*(1-h)-betah*h : 1

    I_K= gK*n**4*(EK-vm) : amp
    alphan = .032*(52*mV+vm)/(1-exp(-(vm+52*mV)/(5*mV)))/mV/ms : Hz
    betan = .5*exp(-(vm+57*mV)/(40*mV))/ms : Hz
    dn/dt = alphan*(1-n)-betan*n : 1

    dvm/dt = membrane_Im/C : volt

    I_Ca = gCa*sinf*(ECa-vm) : amp
    sinf=1/(1+exp(-(vm+25*mV)/(5*mV))) : 1

    I_M = gM*w*(EM-vm) : amp
    dw/dt = (winf-w)/tauw: 1
    tauw = 100*ms : second
    winf = 1/(1+exp(-(vm+25*mV)/(5*mV))) :1


    I_AHP = gAHP*q*(EAHP-vm) : amp
    tauCa = 80*ms : second
    dCa/dt = (3*tauCa/ms*I_Ca/amp-(Ca))/tauCa : 1
    alphaq = (0.02*Ca)/ms : Hz
    betaq = 0.001/ms : Hz
    tauq = 1/(alphaq + betaq)  : second
    qinf =  alphaq/(alphaq + betaq) : 1
    dq/dt = (qinf-q)/tauq: 1

    taum=1/(alpham+betam)/ms : 1
    taun=1/(alphan+betan)/ms : 1
    tauh=1/(alphah+betah)/ms : 1

    minf=alpham/(alpham+betam) :1
    ninf=alphan/(alphan+betan) :1
    hinf=alphah/(alphah+betah) :1
    '''
    #  dq/dt = (qinf-q)/tauq: 1
    #    tauq = (0.0338)/(0.00001*Ca+0.001)*ms : second
    #    qinf = (0.0005*Ca)**2 :1
    #
    #    dCa/dt= (3*I_Ca/amp - 0.0167*Ca) /ms: 1

    #    q = Ca/(Ca+0.5) :1
    #    tauCa = 80*ms : second
    #    dCa/dt = (4*tauCa/ms*I_Ca/amp-(Ca))/tauCa : 1


    #    dCa/dt = (3*tauCa/ms*I_Ca/amp-(Ca))/tauCa : 1
    #    tauCa = 80*ms : second
    #    I_AHP = gAHP*q*(EAHP-vm) : amp
    #    alphaq = (0.000011*Ca)/ms : Hz
    #    betaq = 0.001/ms : Hz
    #    tauq = 1/(alphaq + betaq)  : second
    #    qinf =  alphaq/(alphaq + betaq) : 1
    #    dq/dt = (qinf-q)/tauq: 1
    neuron = b2.NeuronGroup(1, eqs, method='exponential_euler')

    # parameter initialization
    neuron.vm = -66 * b2.mV
    neuron.m = 0.0529324852572
    neuron.h = 0.95
    neuron.n = 0.05
    neuron.Ca = 0

    # tracking parameters
    rec = b2.StateMonitor(neuron,
                          ['vm', 'I_e', 'I_K', 'I_Na', 'I_Ca', 'I_M', 'I_AHP', 'm', 'n', 'h', 'minf', 'ninf',
                           'hinf', 'taum', 'tauh', 'taun', 'Ca', 'sinf', 'winf',

                           'q', 'w'], record=True)

    # running the simulation
    b2.run(simtime)
    return rec

def HH_Step(I_tstart=20, I_tend=180, I_amp=7,
                  tend=200, do_plot=True, model='HH'):

    # 1ms sampled step current
    tmp = np.zeros(np.round(tend)) * b2.uamp
    tmp[int(I_tstart):int(I_tend)] = I_amp * b2.uamp
    curr = b2.TimedArray(tmp, dt=1. * b2.ms)
    # Neuron Model

    if model == "HH":
        rec = HH_Neuron(curr, tend * b2.ms)

        if do_plot:
            plot_data_S(
                rec,
                title="Step current",
            )

    if model == "Adapt":
        rec = HH_Neuron_Adapt(curr, tend * b2.ms)

        if do_plot:
            plot_data_Adapt(
                rec,
                title="Step current",
            )
    return rec

