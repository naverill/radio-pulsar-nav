import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys

# define two fitting functions, linear and 2nd order
def linear_func(x, a, b):
    return a*x + b

def parabolic_func(x, a, b, c):
    return a*x**2 + b*x + c

#############################################################################

# start of the main program

# set up the window size for the plot
plt.figure(figsize=(10,10))

for pol in range(2):

    print("pol= ",pol)
    
    # specify which data file to process
#    filename = 'Vela_2p3_Strong_Pol_040524.csv'
#    filename = 'Vela_2p3_Strong_Pol_140524.csv'
#    filename = 'Vela_2p3_Strong_Pol_010624.csv'
#    filename = "../../Desktop/Vela_2p3_Strong_Pol_120624.csv"

    if pol == 0:
#        filename = "Vela_2p3_Strong_Pol_020724.csv"
        filename = "Vela_2p3_Strong_Pol_171024.csv"
        polename = 'Strong'
        colour = 'b'
        offset = 0.0
    else:
#        filename = "Vela_2p3_Weak_Pol_020724.csv"
        filename = "Vela_2p3_Weak_Pol_171024.csv"
        polename = 'Weak'
        colour = 'g'
        offset = -20.0E-6
        
    # number of microseconds per day
    usecperday = 1E6*60*60*24

    # load the data
    data = np.loadtxt(filename, usecols=(2,3,4,5,8), delimiter=',', dtype=str, skiprows=2)

    # set up empty variables for the analysis
    mjd = []
    snr = []
    p0 = []
    p0_err = []
    df_f = []

    # go through the data line by line, saving valid entries only
    # rows with null values for relevant data are skipped
    # data are saved as strings initially, converted to floats later
    for i in range(len(data)):
        if data[i][0] != '' and data[i][1] != '' and data[i][4] != '':
            mjd.append(data[i][0])
            snr.append(data[i][1])
            p0.append(data[i][2])
            p0_err.append(data[i][3])
            df_f.append(data[i][4])

    # convert all data arrays to floating point        
    mjd = np.array(mjd,dtype=float)
    snr = np.array(snr,dtype=float)
    p0 = np.array(p0,dtype=float)
    p0_err = np.array(p0_err,dtype=float)
    df_f = np.array(df_f,dtype=float)

    # define start MJD -- all data prior to this is dropped
    MJD_start = 60310 

    # define MJD of the glitch
    MJD_glitch = 60430

    # retain only data after the start MJD
    mask = mjd > MJD_start
    mjd = mjd[mask]
    snr = snr[mask]
    p0 = p0[mask]
    p0_err = p0_err[mask]
    df_f = df_f[mask]

    p0 += offset
    
    # recompute the df/f values rather than use CSV values
    f = 1.0/p0
    df = np.zeros(len(f))
    df[1:-1] = f[1:-1] - f[0:-2] 
    df_f = df/f

    for i in range(len(df_f)):
        print(i,df_f[i])
        
    #plt.plot(df_f)
    #plt.show()

    ax1 = plt.subplot(311)
    
    if pol == 0 :
        pref = p0[0]
        
    #plt.plot(mjd, (p0-pref)*1000, colour+'.',label=polename+" pole")
    
    #plt.ylabel("$\Delta P0$ [$\mu$sec]")
    #plt.xlabel("MJD [days]")

    #plt.axvline(MJD_glitch,alpha=0.5,label="Glitch epoch")
    #plt.legend()
    #plt.show()
    #sys.exit()


    # measure how many data points we have after the glitch
    Npostglitch = int(sum(mjd>MJD_glitch))
    print("Number of data points post-glitch = ",Npostglitch)


    #####################################################################
    # first sub plot -- this is SNR versus MJD
    plt.plot(mjd, snr, colour+'.', ms=6,label=polename)
    plt.ylabel("S/N")
    plt.ylim(0,)
    plt.axvline(MJD_glitch,alpha=0.5)
    plt.legend()


    #####################################################################
    # second sub plot -- change in the period versus MJD
    # change in period is measured in usec, relative to the first obs
    ax2 = plt.subplot(312)
    print("Period for first observation = ",p0[0]," ms")
    plt.plot(mjd, (p0-pref)*1000, colour+'.', label=polename)
    plt.ylabel("$\Delta P0$ [$\mu$sec]")
    #plt.xlabel("MJD [days]")
    plt.axvline(MJD_glitch,alpha=0.5)

    # set up xdata and ydata for two fits -- linear and 2nd order
    # pre-glitch data first
    xdata = mjd
    ydata = (p0-pref)*1000
    xdata = xdata[0:-Npostglitch]
    ydata = ydata[0:-Npostglitch]

    # mark the period in the observation prior to the glitch as a horizontal line
    # the idea is to show when the period post-glitch recovers to this value
    #plt.axhline(ydata[-1],alpha=0.2,color='k')

    # do the fit
    popt, pcov = curve_fit(linear_func, xdata, ydata)

    # measure the scatter around the fit
    pre_glitch_scatter = np.std(linear_func(xdata, popt[0], popt[1])-ydata)
    print("Pre-glitch scatter = ",pre_glitch_scatter, " usec")

    # set up some nice limits to show the fit
    # the fit is extended somewhat beyond the last point for
    # comparison with the post-glitch data
    if pol == 0:
        xlo = min(xdata)
        xhi = max(xdata)
        xrange = xhi - xlo
        xlo -= xrange*0.15
        xhi += xrange*1.70

    # draw the pre-glitch 1st order fit
    xvalues = np.linspace(xlo,xhi,100)
    yvalues = linear_func(xvalues, popt[0], popt[1])
    #plt.plot(xvalues, yvalues, 'g-',label=polename + " pole: Pre-glitch")
    print("Pre-glitch Pdot = ",popt[0]," usec/day")
    print("Pre-glitch Pdot = ",popt[0]/usecperday)

    # now set up the post-glitch fit -- linear only
    xdata = mjd
    ydata = (p0-pref)*1000
    xdata = xdata[-Npostglitch+2:]
    ydata = ydata[-Npostglitch+2:]

    popt, pcov = curve_fit(linear_func, xdata, ydata)
    post_glitch_scatter = np.std(linear_func(xdata, popt[0], popt[1])-ydata)
    print("Post-glitch scatter = ",post_glitch_scatter, " usec")
    print("Post-glitch Pdot = ",popt[0]," usec/day")
    print("Post-glitch Pdot = ",popt[0]/usecperday)

    xvalues = np.linspace(xdata[0],xdata[-1],100)
    yvalues = linear_func(xvalues, popt[0], popt[1])

    # draw the post-glitch fit
    #plt.plot(xvalues, yvalues, 'r-',label=polename + ' pole: Post-glitch')

    plt.legend()


    #####################################################################
    # third sub plot -- change in the df/f versus MJD
    ax3 = plt.subplot(313)
    plt.plot(mjd, df_f, colour+'.-',label=polename)
    plt.xlabel("MJD [days]")
    plt.ylabel("$\Delta(\delta f/f)$")
    plt.axvline(MJD_glitch,alpha=0.5)
    plt.legend()


# share x-axes on all plots -- zooming in on any plot affects
# the other two plots in the same way on the x-axis
ax1.sharex(ax2)
ax2.sharex(ax3)
ax3.sharex(ax1)

# set the x-axis limits on the subplots to be the same for all three
plt.subplot(311)
plt.xlim(xlo,xhi)
plt.subplot(312)
plt.xlim(xlo,xhi)
plt.subplot(313)
plt.xlim(xlo,xhi)

plt.show()
plt.savefig("../outputs/vela_snr_v_mjd.png")


print(len(df_f), Npostglitch)

pre_glitch_data = df_f[0:-Npostglitch-1]
post_glitch_data = df_f[-Npostglitch:-1]

print(len(pre_glitch_data))
print(len(post_glitch_data))

return_vals = plt.hist(post_glitch_data,bins=100,label="Post-glitch, N ="+str(len(post_glitch_data)))
plt.hist(pre_glitch_data,bins=return_vals[1],label="Pre-glitch, N ="+str(len(pre_glitch_data)))

#return_vals = plt.hist(pre_glitch_data,bins=100,label="Pre-glitch, N ="+str(len(pre_glitch_data)))
#plt.hist(post_glitch_data,bins=return_vals[1],label="Post-glitch, N ="+str(len(post_glitch_data)))

plt.xlabel("$\Delta(\delta f/f)$")
plt.ylabel("N")

plt.legend()

# plt.show()
plt.savefig("../outputs/vela_post_glitch.png")

