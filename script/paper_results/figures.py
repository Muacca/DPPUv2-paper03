"""Render supplied data without invoking any geometry computation."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .data import load_data

def render(source,output):
    data,scalars=load_data(source)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                         "figure.dpi":150,"savefig.dpi":190,"svg.hashsalt":"paper03-v3"})
    colors=["#17466e","#c46a23","#148576"]
    paths=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ["png","pdf"]:
            path=out/f"{name}.{ext}"
            fig.savefig(path,bbox_inches="tight");paths.append(path)
        plt.close(fig)
    d=data["radius.csv"];fig,ax=plt.subplots(figsize=(7,3.7))
    for name,c in zip(["MX","AX","VT"],colors):ax.plot(d["r"],d[name],label=name,color=c)
    ax.axvline(3,color=".6",ls=":",lw=1);ax.axhline(0,color=".6",lw=.7)
    ax.set(xlabel="Radius scale r (physical round radius r/2)",ylabel="Restricted potential U")
    ax.legend(title="Geometric torsion");save(fig,"fig01_radius")
    fig,axs=plt.subplots(1,2,figsize=(8,3.5))
    for ax,name in zip(axs,["MX","AX"]):
        H=np.array(scalars[name]["hessian"])/np.pi**2
        im=ax.imshow(H,cmap="RdBu_r",vmin=-5,vmax=5)
        ax.set_xticks(range(4),["r","L","eta","V"]);ax.set_yticks(range(4),["r","L","eta","V"])
        ax.set_title(f"{name}: H / pi^2")
        for i in range(4):
            for j in range(4):ax.text(j,i,f"{H[i,j]:g}",ha="center",va="center",color="black")
    save(fig,"fig02_scalar")
    d=data["twist.csv"];fig,ax=plt.subplots(figsize=(7,3.5))
    ax.plot(d["omega1"],d["U"],label="Original frame",color=colors[0])
    ax.plot(d["omega1"],d["U_rotated_frame"],"--",label="Rotated frame",color=colors[1])
    ax.set(xlabel="omega1 (omega2=omega3=0)",ylabel="Restricted potential U")
    ax.legend();save(fig,"fig03_twist")
    d=data["twist_kinetic.csv"];fig,axs=plt.subplots(1,2,figsize=(8,3.5))
    c=np.pi**2/3
    axs[0].plot(d["rho"],d["static"]/c,label="Static: rho^2",color=colors[0])
    axs[0].plot(d["rho"],d["rotating_derivative"]/c,label="Rotating direction: pi^2 rho^4 / 4",color=colors[2])
    axs[0].set(xlabel="Twist amplitude rho",ylabel="Action density / C_omega")
    axs[0].legend(fontsize=8)
    axs[1].plot(d["rho"],d["G_transverse"]/c,label="Transverse (double)",color=colors[0])
    axs[1].plot(d["rho"],np.zeros_like(d["rho"]),ls="--",label="Radial",color=colors[1])
    axs[1].set(xlabel="Twist amplitude rho",ylabel="Velocity Hessian eigenvalue / C_omega")
    axs[1].legend(fontsize=8)
    save(fig,"fig04_kinetic")
    d=data["shear.csv"];fig,axs=plt.subplots(1,2,figsize=(8.5,3.6))
    for name,label,c in zip(["H0","H2","H1"],["m=0","m=+/-2","m=+/-1"],colors):
        axs[0].plot(d["epsilon"],d[name]/(24*np.pi**2),label=label,color=c)
    axs[0].set(xlabel="epsilon",ylabel="Second derivative / (24 pi^2)",ylim=(-3,6),title="Exponential deformation")
    axs[0].legend()
    axs[1].plot(d["epsilon"],d["H0"]/(24*np.pi**2),color=colors[0],label="Exponential")
    axs[1].plot(d["epsilon"],d["H0_linear"]/(24*np.pi**2),color=colors[1],ls="--",label="Normalized linear")
    axs[1].set(xlabel="epsilon",ylabel="Singlet curvature / (24 pi^2)",ylim=(-1,3),title="Dependence on the path")
    axs[1].legend()
    for ax in axs:ax.axhline(0,color=".5",lw=.7);ax.axvline(0,color=".7",lw=.7)
    save(fig,"fig05_shear")
    d=data["bifurcation.csv"];fig,axs=plt.subplots(1,2,figsize=(8,3.4))
    axs[0].plot(d["epsilon"],d["q2_positive"],color=colors[0])
    axs[0].set(xlabel="epsilon (fixed control)",ylabel="Doublet amplitude rho at the minimum")
    axs[1].plot(d["epsilon"],d["delta_U"],color=colors[2]);axs[1].set(xlabel="epsilon",ylabel="U(branch) - U(q2=0)")
    for ax in axs:ax.axvline(1/np.sqrt(2)-1,color=".5",ls=":",lw=1)
    save(fig,"fig06_branch")
    d=data["shape_minima.csv"];fig,ax=plt.subplots(figsize=(7,3.6))
    ax.plot(d["epsilon"],d["F_axis"],"--",color=".6",label="Axisymmetric")
    ax.plot(d["epsilon"],d["F_min"],color=colors[0],label="Doublet minimized")
    ax.plot(d["epsilon"],d["F_global"],color=colors[1],label="Full shape minimum")
    ax.axvline(1/np.sqrt(2)-1,color=".5",ls=":",lw=1)
    ax.set(xlabel="Axial control epsilon",ylabel="Dimensionless shape potential F",ylim=(-3.2,.5))
    ax.legend(fontsize=9);save(fig,"fig07_shape_minima")
    d=data["invariants.csv"];fig,ax=plt.subplots(figsize=(7,3.5))
    ax.plot(d["q1"],d["delta_U"],label="Exact potential",color=colors[0])
    ax.plot(d["q1"],d["quartic"],"--",label="Through fourth order",color=colors[2])
    ax.set(xlabel="q1 on h=q1 T1",ylabel="U(h) - U(0)");ax.legend()
    save(fig,"fig08_expansion")
    d=data["cross_kinetic.csv"];cdata=data["coupled_action.csv"]
    fig,axs=plt.subplots(1,2,figsize=(9.2,3.7))
    for name,label,c in zip(["G_rw","G_sw","G_rs"],["G(r,w)","G(s,w)","G(r,s)"],colors):
        axs[0].plot(d["w"],d[name]/np.pi**2,color=c,label=label)
    axs[0].axhline(0,color=".6",lw=.7)
    axs[0].set(xlabel="One-axis twist w",ylabel="Velocity cross coefficient / pi^2",title="Derivative couplings")
    axs[0].legend(fontsize=8)
    axs[1].plot(cdata["n"],cdata["delta_static"],"--",color=colors[1],label="Static contribution")
    axs[1].plot(cdata["n"],cdata["derivative"],":",color=colors[2],label="Derivative contribution")
    axs[1].plot(cdata["n"],cdata["delta_total"],"o-",markersize=3,color=colors[0],label="Total relative to round")
    axs[1].axhline(0,color=".6",lw=.7)
    axs[1].set(xlabel="Circle frequency n",ylabel="Action contribution",title="Coupled periodic family")
    axs[1].legend(fontsize=8)
    save(fig,"fig09_coupled")
    d=data["curvature.csv"];fig,axs=plt.subplots(1,2,figsize=(8,3.4))
    axs[0].plot(d["r"],d["B"],color=colors[0],label="Internal pseudoscalar B")
    axs[0].axhline(0,color=".7",lw=.7);axs[0].set(xlabel="r",ylabel="B",title="Internal dual")
    axs[1].plot(d["r"],d["P"],color=colors[2]);axs[1].set(xlabel="r",ylabel="Standard P",title="Form-index contraction",ylim=(-1e-12,1e-12))
    axs[1].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    save(fig,"fig10_curvature")
    return paths
