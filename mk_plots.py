import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import romberg
import h5py
from scipy.interpolate import interp1d
import sys
import os

delta_vir=200.
rhocrit0=27.755
lbox=210
nside=3072
nbins=100
lmmin=7
lmmax=13.3

#### Define functions

# Compute the quantity E(a) in H(a)=H0E(a)
def eofa(aexp,omegam,omegal):
    return np.sqrt(omegam/aexp**3+omegal)
    
# Read snapshot information from the .siminfo file
def read_infofile(infile):
    descriptor,value=np.loadtxt(infile,dtype={'names':('descriptor','value'),'formats':('S20','f')},delimiter=':',unpack=True)

    for i in range(len(descriptor)):
        if 'h_val' in descriptor[i].decode():
            hubble=value[i]
        if 'ScaleFactor' in descriptor[i].decode():
            expansion=value[i]
        if 'Omega_m' in descriptor[i].decode():
            omegam=value[i]
        if 'Omega_Lambda' in descriptor[i].decode():
            omegal=value[i]

    return omegam,omegal,hubble,expansion

# Read halo bulk properties from the .properties file
def read_properties(dir_path,filename_root,verbose):
    dir_list=os.listdir(dir_path)
    files=[]
    for i in range(len(dir_list)):
        if filename_root in dir_list[i]:
            files+=[dir_list[i]]
    print("Reading "+str(len(files))+" files in "+dir_path )

    m200=[]
    mfof=[]
    rmax=[]
    r200=[]
    vmax=[]
    spin=[]
    nsub=[]
    hid=[]
    for i in range(len(files)-1):
        if verbose==True:
            print("Reading data from "+files[i])
        f=h5py.File(dir_path+files[i],"r")
        mass_200crit=f['Mass_200crit'][()]
        mass_fof=f['Mass_FOF'][()]
        numSubStruct=f['numSubStruct'][()]
        hostHaloID=f['hostHaloID'][()]
        rvcmax=f['Rmax'][()]
        vcmax=f['Vmax'][()]
        lambdab=f['lambda_B'][()]
        radius_200crit=f['R_200crit'][()]
        f.close()
        # Append values to array
        m200=np.append(m200,mass_200crit)
        mfof=np.append(mfof,mass_fof)
        nsub=np.append(nsub,numSubStruct)
        hid=np.append(hid,hostHaloID)
        r200=np.append(r200,radius_200crit)
        rmax=np.append(rmax,rvcmax)
        vmax=np.append(vmax,vcmax)
        spin=np.append(spin,lambdab)
    return m200,mfof,nsub,hid,r200,rmax,vmax,spin

def get_mass_function(mass,nbins,logmass_min,logmass_max,lbox):
    dlm=(lmmax-lmmin)/nbins
    (n,lm)=np.histogram(np.log10(mass)+10.,bins=nbins,range=([logmass_min,logmass_max]))
    logmass_bins=np.zeros(len(n))
    for i in range(len(lm)-1):
        logmass_bins[i]=np.log10(0.5*(10**lm[i]+10**lm[i+1]))
    dndlm=n/dlm/lbox**3
    return logmass_bins,dndlm

def setup_plot(left,bottom,width,height,spacing,xmin,xmax,ymin,ymax):
    # definitions for the axes
    rect_upper = [left, bottom + 0.2 + spacing, width, height]
    rect_lower = [left, bottom                , width, 0.2]

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_axes(rect_upper)
    ax.tick_params(axis="y", labelleft=True)
    ax.tick_params(axis="x", labelbottom=False)
    ax.grid(True)
    plt.ylabel('Log$_{10}$ dn/dLog$_{10}$M [$h^3$ Mpc$^{-3}$]')
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])

    ax_lower = fig.add_axes(rect_lower, sharex=ax)
    ax_lower.grid(True)
    plt.ylim([-1,1])
    plt.ylabel('$\Delta$ Log$_{10}$ dn/dLog$_{10}$M')
    plt.xlabel('Log$_{10}$ M$_{200}$ [$h^{-1}$ M$_{\odot}$]')
    
    return ax,ax_lower

def plot_single_mass_function(dir_path,snap_num,nside):
    filename="snap_%04d.VELOCIraptor.siminfo"%int(snap_num)

    omegam,omegal,HubbleParam,aexp=read_infofile(dir_path+filename)

    print("HubbleParam: "+str(HubbleParam))
    print("Omega0: "+str(omegam))
    print("OmegaLambda: "+str(omegal))
    print("Redshift: "+str(1./aexp-1))
    
    lmp=np.log10(omegam*rhocrit0*(lbox/nside)**3)+10
    
    print("Particle mass [Log10 Solar masses]: ",(lmp-np.log10(HubbleParam)))

    filename_root="snap_%04d.VELOCIraptor.properties"%int(snap_num)

    m200,mfof,nsub,hid=read_properties(dir_path,filename_root,verbose=True)

    # Apply cuts to the data
    mvir=HubbleParam*m200[(hid==-1)&(m200>0)]
    mtot=HubbleParam*mfof[(hid==-1)&(mfof>0)]
    msub=HubbleParam*m200[(hid!=-1)&(m200>0)]

    # Get the binned mass function data
    logmass_bins,dndlm=get_mass_function(mvir,nbins,lmmin,lmmax,lbox)

    # Read HMFcalc data and interpolate
    (m_b,n_b)=np.loadtxt("./behroozi_z5_default.txt",usecols=(0,3),unpack=True)
    (m_w,n_w)=np.loadtxt("./watson_z5_default.txt",usecols=(0,3),unpack=True)

    f_int_b=interp1d(np.log10(m_b),np.log10(n_b),kind="cubic")
    f_int_w=interp1d(np.log10(m_w),np.log10(n_w),kind="cubic")
    f_b=f_int_b(logmass_bins)
    f_w=f_int_w(logmass_bins)

    # Set up plot
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.0
    xmin=7
    xmax=13.3
    ymin=-6.5
    ymax=1

    ax1,ax2=setup_plot(left,bottom,width,height,spacing,xmin,xmax,ymin,ymax)

    # Now add data to the plots
    ax1.plot(logmass_bins,np.log10(dndlm),'k',color='black',lw=3)
    ax1.plot(np.log10(m_b),np.log10(n_b),'k:',color='red',lw=3)
    ax1.plot(np.log10(m_w),np.log10(n_w),'k--',color='blue',lw=3)
    ax1.legend(('VELOCIraptor M$_{200}$','Behroozi et al. 2013','Watson et al. 2013'),loc="upper right")
    ax1.axvline(lmp+1,color='black',ls=':')
    ax1.axvline(lmp+2,color='black',ls=':')
    ax1.axvline(lmp+3,color='black',ls=':')
    ax1.text(lmp+0.8,-6,'10 particles',rotation='vertical')
    ax1.text(lmp+1.8,-6,'100 particles',rotation='vertical')
    ax1.text(lmp+2.8,-6,'1000 particles',rotation='vertical')

    ax2.plot(logmass_bins,np.log10(dndlm)-f_b,'k:',color="red",lw=3)
    ax2.plot(logmass_bins,np.log10(dndlm)-f_w,'k--',color="blue",lw=3)
    ax2.axhline(np.log10(0.5),color='black',ls='--',lw=1)
    ax2.axhline(np.log10(1.5),color='black',ls='--',lw=1)
    ax2.text(8,np.log10(0.5)-0.15,'-50%',rotation='horizontal')
    ax2.text(8,np.log10(1.5)+0.1,'+50%')

    # Finished adding data to the plot, so write to file
    filename="snap_%04d.SO200_mass_function.pdf"%int(snap_num)
    plt.savefig(dir_path+filename, bbox_inches='tight')

def plot_multiple_mass_functions(file_list,nside_arr,lbox_arr,color_arr,lw_arr,box_params,mf_params,outfile,verbose):
    print(file_list)
    # First set up the plotting area
    left, width, bottom, height, spacing = box_params
    xmin, xmax, ymin, ymax = mf_params
    print(box_params)
    print(mf_params)
    ax1,ax2=setup_plot(left,bottom,width,height,spacing,xmin,xmax,ymin,ymax)
 
    # Read HMFcalc data
    (m_b,n_b)=np.loadtxt("/Users/cpower/Downloads/behroozi_z5_default.txt",usecols=(0,3),unpack=True)
    (m_w,n_w)=np.loadtxt("/Users/cpower/Downloads/watson_z5_default.txt",usecols=(0,3),unpack=True)

    # Now read in the data...
    om=np.zeros(len(file_list))
    ol=np.zeros(len(file_list))
    hub=np.zeros(len(file_list))
    aexp=np.zeros(len(file_list))

    for i in range(len(file_list)):
        filename=file_list[i]+'.siminfo'
        om[i],ol[i],hub[i],aexp[i]=read_infofile(filename)
        if verbose==True:
            print("Filename: ",filename)
            print("HubbleParam: "+str(hub[i]))
            print("Omega0: "+str(om[i]))
            print("OmegaLambda: "+str(ol[i]))
            print("Redshift: "+str(1./aexp[i]-1))

        lmp=np.log10(om[i]*rhocrit0*(lbox_arr[i]/nside_arr[i])**3)+10
    
        if verbose==True:
            print("Particle mass [Log10 Solar masses]: ",(lmp-np.log10(hub[i])))

        dir_path,filename_root=os.path.split(file_list[i])
        dir_path=dir_path+'/'
        filename_root=filename_root+'.properties'

        m200,mfof,nsub,hid=read_properties(dir_path,filename_root,verbose=True)

        # Apply cuts to the data
        mvir=hub[i]*m200[(hid==-1)&(m200>0)]
        mtot=hub[i]*mfof[(hid==-1)&(mfof>0)]
        msub=hub[i]*m200[(hid!=-1)&(m200>0)]

        # Get the binned mass function data
        logmass_bins,dndlm=get_mass_function(mvir,nbins,lmmin,lmmax,lbox)
 
        # Interpolate the mass functions to the theoretical comparisons
        f_int_b=interp1d(np.log10(m_b),np.log10(n_b),kind="cubic")
        f_int_w=interp1d(np.log10(m_w),np.log10(n_w),kind="cubic")
        f_b=f_int_b(logmass_bins)
        f_w=f_int_w(logmass_bins)
 
        # Now add data to the plots
        ax1.plot(logmass_bins,np.log10(dndlm),'k',color=color_arr[i],lw=lw_arr[i])
        ax1.plot(np.log10(m_b),np.log10(n_b),'k:',color='black',lw=2)
        ax1.plot(np.log10(m_w),np.log10(n_w),'k--',color='black',lw=2)
        ax2.plot(logmass_bins,np.log10(dndlm)-f_b,'k:',color=color_arr[i],lw=lw_arr[i])
        ax2.plot(logmass_bins,np.log10(dndlm)-f_w,'k--',color=color_arr[i],lw=lw_arr[i])

        ax1.axvline(lmp+1.8,color='black',ls=':')
        ax1.text(lmp+1.6,-6,'~60 particles',rotation='vertical')
 
        ax1.axvline(lmp+3,color='black',ls=':')
        ax1.text(lmp+2.8,-6,'1000 particles',rotation='vertical')
        ax2.axvline(lmp+3,color='black',ls=':')

        if i==0:
            ax1.legend(('_','Behroozi et al. 2013','Watson et al. 2013'),loc="upper right")

    ax2.axhline(np.log10(0.5),color='black',ls='--',lw=1)
    ax2.axhline(np.log10(1.5),color='black',ls='--',lw=1)
    ax2.text(8,np.log10(0.5)-0.15,'-50%',rotation='horizontal')
    ax2.text(8,np.log10(1.5)+0.1,'+50%')
    ax1.text(12.5,-0.7,'z=5',color='black',weight='heavy',size='large')

    # Finished adding data to the plot, so write to file
    plt.savefig(outfile, bbox_inches='tight')

def plot_multiple_rmax_vmax(file_list,nside_arr,lbox_arr,color_arr,lw_arr,box_params,mf_params,outfile,verbose):

    # First set up the plotting area
    left, width, bottom, height, spacing = box_params
    xmin, xmax, ymin, ymax = mf_params

    ax1,ax2=setup_plot(left,bottom,width,height,spacing,xmin,xmax,ymin,ymax)
 
    # Now read in the data...
    om=np.zeros(len(file_list))
    ol=np.zeros(len(file_list))
    hub=np.zeros(len(file_list))
    aexp=np.zeros(len(file_list))

    for i in range(len(file_list)):
        filename=file_list[i]+'.siminfo'
        om[i],ol[i],hub[i],aexp[i]=read_infofile(filename)
        if verbose==True:
            print("Filename: ",filename)
            print("HubbleParam: "+str(hub[i]))
            print("Omega0: "+str(om[i]))
            print("OmegaLambda: "+str(ol[i]))
            print("Redshift: "+str(1./aexp[i]-1))

        lmp=np.log10(om[i]*rhocrit0*(lbox_arr[i]/nside_arr[i])**3)+10
    
        if verbose==True:
            print("Particle mass [Log10 Solar masses]: ",(lmp-np.log10(hub[i])))

        dir_path,filename_root=os.path.split(file_list[i])
        dir_path=dir_path+'/'
        filename_root=filename_root+'.properties'

        m200,mfof,nsub,hid,r200,rvcmax,vcmax,lmbda=read_properties(dir_path,filename_root,verbose=True)

        # Apply cuts to the data
        rmax=hub[i]*rvcmax[(hid==-1)&(m200>0)]/aexp
        vmax=hub[i]*vcmax[(hid==-1)&(m200>0)]

        # Get the binned mass function data
        logmass_bins,dndlm=get_mass_function(mvir,nbins,lmmin,lmmax,lbox)
 
        # Now add data to the plots
        ax1.plot(logmass_bins,np.log10(dndlm),'k',color=color_arr[i],lw=lw_arr[i])
        ax1.plot(np.log10(m_b),np.log10(n_b),'k:',color='black',lw=2)
        ax1.plot(np.log10(m_w),np.log10(n_w),'k--',color='black',lw=2)
        ax2.plot(logmass_bins,np.log10(dndlm)-f_b,'k:',color=color_arr[i],lw=lw_arr[i])
        ax2.plot(logmass_bins,np.log10(dndlm)-f_w,'k--',color=color_arr[i],lw=lw_arr[i])

        ax1.axvline(lmp+1.8,color='black',ls=':')
        ax1.text(lmp+1.6,-6,'~60 particles',rotation='vertical')
 
        ax1.axvline(lmp+3,color='black',ls=':')
        ax1.text(lmp+2.8,-6,'1000 particles',rotation='vertical')
        ax2.axvline(lmp+3,color='black',ls=':')

        if i==0:
            ax1.legend(('_','Behroozi et al. 2013','Watson et al. 2013'),loc="upper right")

    ax2.axhline(np.log10(0.5),color='black',ls='--',lw=1)
    ax2.axhline(np.log10(1.5),color='black',ls='--',lw=1)
    ax2.text(8,np.log10(0.5)-0.15,'-50%',rotation='horizontal')
    ax2.text(8,np.log10(1.5)+0.1,'+50%')
    ax1.text(12.5,-0.7,'z=5',color='black',weight='heavy',size='large')

    # Finished adding data to the plot, so write to file
    plt.savefig(dir_path+outfile, bbox_inches='tight')

####



def main(argv):
    dir_path=argv[0]
    snap_num=argv[1]

    nside=3072
    plot_single_mass_function(dir_path,snap_num,nside)

#    nside_arr=[1536,3072]
#    lbox_arr=[210,210]
#    color_arr=['blue','red']
#    lw_arr=[3,3]
#    file_list=(dir_path+'L210N1536/snap_0119.VELOCIraptor',dir_path+'L210N3072/snap_0119.VELOCIraptor')
#    box_params=[0.1,0.65,0.1,0.65,0]  ## left,width,bottom,height,spacing
#    mf_params=[8.7,13.3,-6.3,0.3]     ## lmmin,lmmax,ldndlm_min,ldndlm_max

#    filename="mass_function_l210_n1536_n3072_z5.pdf"

#    plot_multiple_mass_functions(file_list,nside_arr,lbox_arr,color_arr,lw_arr,box_params,mf_params,filename,verbose=True)

#   mf_params=[8.7,13.3,-6.3,0.3]     ## lmmin,lmmax,ldndlm_min,ldndlm_max

#    filename="vmax_rmax_l210_n1536_n3072_z5.pdf"

#    plot_multiple_rmax_vmax(file_list,nside_arr,lbox_arr,color_arr,lw_arr,box_params,mf_params,filename,verbose=True)

#    mf_params=[8.7,13.3,-6.3,0.3]     ## lmmin,lmmax,ldndlm_min,ldndlm_max

#    filename="mass_function_l210_n1536_n3072_z5.pdf"

#    plot_multiple_spins(file_list,nside_arr,lbox_arr,color_arr,lw_arr,box_params,mf_params,filename,verbose=True)

if __name__ == "__main__":
    print(sys.argv[1:])
    main(sys.argv[1:])

