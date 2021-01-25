
import fileinput
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from matplotlib import cm
import numpy as np
import glob
import os


PATH = r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\G1_grand_selected"
CPATH = r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\G1_constraint_100"
PATH2 = r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data"
PATHrif = r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\rif"
chromosome1_f_e = pd.DataFrame(np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr1_f_e.npy"))
Chromosome1_r_d = pd.DataFrame(np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr1_r_d.npy"))
Chromosome2_f_e = np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr2_f_e.npy")
Chromosome2_r_d = np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr2_r_d.npy")
Chromosome3_f_e = np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr3_f_e.npy")
Chromosome3_r_d = np.load(r"C:\Users\olive\OneDrive\Documents\University\Dissertation\data\Chromosome data\chr3_r_d.npy")



def timing1():
    ef_norm,dr_norm = norm_mean(chromosome1_f_e), norm_mean(Chromosome1_r_d)
    t1 = (ef_norm + dr_norm) / 2.

    return t1

def timing2():
    ef_norm,dr_norm = norm_mean(Chromosome2_f_e), norm_mean(Chromosome2_r_d)
    t2 = (ef_norm + dr_norm) / 2.

    return t2

def timing3():
    ef_norm,dr_norm = norm_mean(Chromosome3_f_e), norm_mean(Chromosome3_r_d)
    t3 = (ef_norm + dr_norm) / 2.

    return t3

def plot_timing():

    plt.figure()
    fig,(ax1, ax2, ax3) = plt.subplots(3,1)
    fig.suptitle("Replication timing in S. pombe")
    ax1.plot(timing1())
    ax2.plot(timing2())
    ax3.plot(timing3())
    ax1.axis([0, 3100, 0, 1])
    ax2.axis([0, 3100, 0, 1])
    ax3.axis([0, 3100, 0, 1])
    ax1.set(xlabel="Position on Chromosome 1 (kb)")
    ax2.set(xlabel="Position on Chromosome 2 (kb)",ylabel = "Replication Time")
    ax3.set(xlabel="Position on Chromosome 3 (kb)")
    plt.tight_layout()

def get_orieff():
    
    return [ np.load(PATH2+"\chr"+str(i)+"_orieff.npy") for i in range(1,4) ]


def plot_orieff2d(source):
    """
    Plots the 2D projection of timing using the files from source (selected or constraint), and the index of the file idx (note: idx is a list of 2 strings for "selected"
    """
    plt.figure()
    if source == "selected":
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))
            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_orieff = get_orieff()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_orieff,lengths,chs_position) ])

            chs_cm = [ cm.hot(a) for a in chs_remap ]
            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in [ch1,ch2,ch3] ]
            for pos,col in zip(positions_2D,chs_cm):
                    plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha = 0.3)
            
            
        plt.colorbar(cm.ScalarMappable(cmap=cm.hot))
        #Nuclear membrane
        T = np.linspace(0,2*np.pi,200)
        X,Y = 1.33*np.cos(T),1.33*np.sin(T)

        #Nucleolus
        alpha = np.arccos(1.51/2 / 1.33 )
        Tn = np.linspace( np.pi - alpha, np.pi + alpha, 200 )
        Xn,Yn = 1.51+1.33*np.cos(Tn), 1.33*np.sin(Tn)

        #SPB
        Xs,Ys = -1.13 + 0.2*np.cos(T), 0.2*np.sin(T)

        #Plot
        plt.plot(X,Y,c="k")
        plt.plot(Xn,Yn,c="k")
        plt.plot(Xs,Ys,c="k")
        plt.show(block=False)
            
            
        

        plt.show(block=False)
    else:
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(cload_data(filename))
            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_orieff = get_orieff()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_orieff,lengths,chs_position) ])
            chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in chs_position ]
            for pos,col in zip(positions_2D,chs_cm):
                    plt.scatter(pos[:,0],pos[:,1],marker=".",c=col)
            
            
        plt.colorbar(cm.ScalarMappable(cmap=cm.hot))

        plt.show(block=False)   
        

def norm_mean(a):
    return 0.5*a / np.nanmean(a)

def timing(ef,dr):
    ef_norm,dr_norm = norm_mean(ef), norm_mean(dr)
    avg = (ef_norm + dr_norm) / 2.
    avg = 2*avg - 1
    avg = 0.5 * avg / np.nanmean(avg)
    return np.nancumsum(avg)


def read_datarif(chro,strand,pol):
    return np.load(PATHrif+"\\rif1d_"+"chr"+str(chro)+"_"+strand+"_"+pol+".npy")

def read_data(chro,strand,pol):
    return np.load(PATH2+"\\chr"+str(chro)+"_"+strand+"_"+pol+".npy")

def get_timingsrif():
    res = []
    for chro in [1,2,3]:
        ef = read_datarif(chro,"f","e")
        dr = read_datarif(chro,"r","d")
        res.append(timing(ef,dr))
    return res


def get_timings():
    res = []
    for chro in [1,2,3]:
        ef = read_data(chro,"f","e")
        dr = read_data(chro,"r","d")
        res.append(timing(ef,dr))

    return res


def normalise(l):
    """
    given a list of arrays, normalises all of them to [0,1], using the max and min of all the arrays
    """
    maxs = [ np.nanmax(a) for a in l ]
    mins = [ np.nanmin(a) for a in l ]

    m = np.min(mins)
    M = np.max(maxs)

    return [ (a - m) / (M - m) for a in l ] 

def normalise2(l):
    """
    given a list of arrays, normalises all of them to [0,1],
    """
    maxs = [ np.nanmax(a) for a in l ]
    mins = [ np.nanmin(a) for a in l ]

    
    
    return [ (a - m) / (M - m) for a,m,M in zip(l,mins,maxs) ] 



def map_scales(a,N):
    """
    Maps array a (len(a) > N) into N bins.
    The result contains N bins with the average values of the corresponding data in a.
    """
    res, counts = np.zeros(N), np.zeros(N)

    #How many positions from a go into a bin in the result
    bins = len(a) // N  + 1 #Note: +1 for the remainder.
    print(bins)
    for i,x in enumerate(a):
        res[ i // bins ] += a[i]
        counts[ i // bins ] += 1.

    return res / counts

def map_scales2(X,Y,Z,f):
    """
    Takes the values Y, at points X, bins them to the closest points in Z and aggregates them using f.

    Arguments:
        X : np array of positions
        Y : np array of values (must have len(X) == len(Y) )
        Z : np array of new positions
        f : aggregate function. Takes a list and returns a value.

    Output:
        array of length len(Z)
    """

    #Find the position in Z closest to each x.
    #Note: very inefficient. We could impose that X and Z are sorted and speed up this part. 
    bins = np.array([ np.argmin( np.abs ( Z - x ) ) for x in X ]) 

    agg = [ Y[ bins == i ] for i,_ in enumerate(Z) ]

    return np.array( [f(a) for a in agg ] )


def load_data(filename):
    """
    Loads data with the G1_grand_i_j name format
    """
    df = pd.read_csv(filename)
    df.columns = ["ch","x","y","z"]
    return df

def load_datasingle(i,j):
    """
    Loads data with the G1_grand_i_j name format
    """
    df = pd.read_csv(PATH+"\G1_grand_"+i+"_"+j+".csv")
    df.columns = ["ch","x","y","z"]
    return df

def cload_data(filename):
    """
    Loads data with the G1_grand_i_j name format
    """
    df = pd.read_csv(filename)
    df.columns = ["ch","x","y","z"]
    return df


def split_chromosomes(df):
    """
    Splits chromosomes in three separate dataframes
    """
    chs = [ df.loc[ df["ch"] ==c ] for c in [1,2,3] ]

    return chs


def cplot(i):
    """
    Plots file with the G1_grand_i_j name format
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    ch1,ch2,ch3 = split_chromosomes(cload_data(i))

    ax.plot(ch1["x"],ch1["y"],ch1["z"],label = "chr. 1")
    ax.plot(ch2["x"],ch2["y"],ch2["z"],label = "chr. 2")
    ax.plot(ch3["x"],ch3["y"],ch3["z"],label = "chr. 3")
    plt.legend()

    plt.show(block=False)

### Projection to 2D plots

def proj(x):
    """
    Assuming symmetry around x axis, projects a 3D point to a 2D representation.
    Uses sign(y) to have a circle plot (rather than half circle)
    """

    return np.array([ x[0], np.sign(x[1])*np.linalg.norm(x[1:]) ] )





def plot_time2d(source):
    """
    Plots the 2D projection of timing using the files from source (selected or constraint)"
    """
    plt.figure()
    if source == "selected":
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))

    
            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timings()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in [ch1,ch2,ch3] ]
            for pos,col in zip(positions_2D,chs_cm):
               plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha = 0.5)

            
            
        plt.colorbar(cm.ScalarMappable(cmap=cm.hot))
        #Nuclear membrane
        T = np.linspace(0,2*np.pi,200)
        X,Y = 1.33*np.cos(T),1.33*np.sin(T)

        #Nucleolus
        alpha = np.arccos(1.51/2 / 1.33 )
        Tn = np.linspace( np.pi - alpha, np.pi + alpha, 200 )
        Xn,Yn = 1.51+1.33*np.cos(Tn), 1.33*np.sin(Tn)

        #SPB
        Xs,Ys = -1.13 + 0.2*np.cos(T), 0.2*np.sin(T)

        #Plot
        plt.plot(X,Y,c="k")
        plt.plot(Xn,Yn,c="k")
        plt.plot(Xs,Ys,c="k")
        plt.show(block=False)

    else:
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))

    
        

            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timings()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in chs_position ]
            for pos,col in zip(positions_2D,chs_cm):
               plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha = 0.5)
                
        
        #Nuclear membrane
        T = np.linspace(0,2*np.pi,200)
        X,Y = 1.33*np.cos(T),1.33*np.sin(T)

        #Nucleolus
        alpha = np.arccos(1.51/2 / 1.33 )
        Tn = np.linspace( np.pi - alpha, np.pi + alpha, 200 )
        Xn,Yn = 1.51+1.33*np.cos(Tn), 1.33*np.sin(Tn)

        #SPB
        Xs,Ys = -1.13 + 0.2*np.cos(T), 0.2*np.sin(T)

        #Plot
        plt.plot(X,Y,c="k")
        plt.plot(Xn,Yn,c="k")
        plt.plot(Xs,Ys,c="k")
        plt.show(block=False)

def plot_time2dRif1(source):
    """
    Plots the 2D projection of timing using the files from source (selected or constraint), and the index of the file idx (note: idx is a list of 2 strings for "selected"
    """
    plt.figure()
    if source == "selected":
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))

    
        

            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timingsrif()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in chs_position ]
            for pos,col in zip(positions_2D,chs_cm):
               plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha = 0.3)
        T = np.linspace(0,2*np.pi,200)
        X,Y = 1.33*np.cos(T),1.33*np.sin(T)
        plt.colorbar(cm.ScalarMappable(cmap=cm.hot))
        #Nucleolus
        alpha = np.arccos(1.51/2 / 1.33 )
        Tn = np.linspace( np.pi - alpha, np.pi + alpha, 200 )
        Xn,Yn = 1.51+1.33*np.cos(Tn), 1.33*np.sin(Tn)

        #SPB
        Xs,Ys = -1.13 + 0.2*np.cos(T), 0.2*np.sin(T)

        #Plot
        plt.plot(X,Y,c="k")
        plt.plot(Xn,Yn,c="k")
        plt.plot(Xs,Ys,c="k")
        plt.show(block=False)

    else:
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))

    
        

            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timingsrif()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ proj(x) for x in pos[["x","y","z"]].to_numpy() ]) for pos in chs_position ]
            for pos,col in zip(positions_2D,chs_cm):
                plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha = 0.5)
            
            
        plt.colorbar(cm.ScalarMappable(cmap=cm.hot))
                #Nuclear membrane
        T = np.linspace(0,2*np.pi,200)
        X,Y = 1.33*np.cos(T),1.33*np.sin(T)

        #Nucleolus
        alpha = np.arccos(1.51/2 / 1.33 )
        Tn = np.linspace( np.pi - alpha, np.pi + alpha, 200 )
        Xn,Yn = 1.51+1.33*np.cos(Tn), 1.33*np.sin(Tn)

        #SPB
        Xs,Ys = -1.13 + 0.2*np.cos(T), 0.2*np.sin(T)

        #Plot
        plt.plot(X,Y,c="k")
        plt.plot(Xn,Yn,c="k")
        plt.plot(Xs,Ys,c="k")
        plt.show(block=False)

def plot(i,j):
    """
    Plots file with the G1_grand_i_j name format
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    ch1,ch2,ch3 = split_chromosomes(load_data(i,j))

    ax.plot(ch1["x"],ch1["y"],ch1["z"],label = "chr. 1")
    ax.plot(ch2["x"],ch2["y"],ch2["z"],label = "chr. 2")
    ax.plot(ch3["x"],ch3["y"],ch3["z"],label = "chr. 3")
    plt.legend()

    plt.show(block=False)

def plot_time(i,j):
    """
    Plots file with the G1_grand_i_j name format
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    lengths = [ 5579133, 4539804, 2452883 ]

    ch1,ch2,ch3 = split_chromosomes(load_datasingle(i,j))
    chs_position = [ch1,ch2,ch3]
    chs_timings = get_timings()
    chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
    
    chs_cm = [ cm.hot(a) for a in chs_remap ]

    for pos,col in zip(chs_position,chs_cm):
        ax.scatter(pos["x"],pos["y"],pos["z"],marker=".",c=col)

    plt.colorbar(cm.ScalarMappable(cmap=cm.hot))

    plt.show(block=False)

def plot_orieff(i,j):
    """
    Plots file with the G1_grand_i_j name format
    """
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    lengths = [ 5579133, 4539804, 2452883 ]
    ch1,ch2,ch3 = split_chromosomes(load_datasingle(i,j))
    chs_position = [ch1,ch2,ch3] 
    chs_orieff= get_orieff()
    chs_remap = normalise2( [ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_orieff,lengths,chs_position) ] )
    chs_cm = [ cm.hot(a) for a in chs_remap ]

    for pos,col in zip(chs_position,chs_cm):
        ax.scatter(pos["x"],pos["y"],pos["z"],marker=".",c=col)

    plt.colorbar(cm.ScalarMappable(cmap=cm.hot))

    plt.show(block=False)

def plot_time_slice(source):
    """
    Plots the 2D projection of timing using the files from source (selected or constraint), and the index of the file idx (note: idx is a list of 2 strings for "selected"
    """
    plt.figure()
    if source == "selected":
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))
            
            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timings()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            #chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ (x[1],x[2]) for x in pos[["x","y","z"]].to_numpy() if -0.5 < x[0] and x[0] < -0.4]) for pos in chs_position ]

            chs_remap2 = [ np.array( [c for x,c in zip( pos[["x","y","z"]].to_numpy(),colors) if -0.5 < x[0] and x[0] < -0.4] ) for pos,colors in zip(chs_position,chs_remap)]
            chs_cm = [ cm.hot(a) for a in chs_remap2 ]
            for pos,col in zip(positions_2D,chs_cm):
                plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha=0.5)
        
    
    plt.show(block=False)
    
def plot_time_slicerif1(source):
    """
    Plots the 2D projection of timing using the files from source (selected or constraint), and the index of the file idx (note: idx is a list of 2 strings for "selected"
    """
    plt.figure()
    if source == "selected":
        for filename in glob.glob(os.path.join(PATH, "*.csv")):
            ch1,ch2,ch3 = split_chromosomes(load_data(filename))
            
            lengths = [ 5579133, 4539804, 2452883 ]

            chs_position = [ch1,ch2,ch3] 
            chs_timings = get_timingsrif()
            chs_remap = normalise2([ map_scales2(np.arange(150,len(Y)*300,300),Y,np.linspace(0,l,len(pos)),np.mean) for Y,l,pos in zip(chs_timings,lengths,chs_position) ])
            #chs_cm = [ cm.hot(a) for a in chs_remap ]

            positions_2D = [ np.array([ (x[1],x[2]) for x in pos[["x","y","z"]].to_numpy() if -0.5 < x[0] and x[0] < -0.4]) for pos in chs_position ]

            chs_remap2 = [ np.array( [c for x,c in zip( pos[["x","y","z"]].to_numpy(),colors) if -0.5 < x[0] and x[0] < -0.4] ) for pos,colors in zip(chs_position,chs_remap)]
            chs_cm = [ cm.hot(a) for a in chs_remap2 ]
            for pos,col in zip(positions_2D,chs_cm):
                plt.scatter(pos[:,0],pos[:,1],marker=".",c=col,alpha=0.5)

    
    plt.show(block=False)




plot_orieff2d("selected")