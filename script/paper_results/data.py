"""Explicit computation of the paper's fixed datasets."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from .model import cartan,potential,shape_matrix,shear_curvatures,shear_branch,scalar_data,shear_basis,exp_symmetric,shape_minimum,twist_bulk
from .circle import diagonal_cross_coefficients,coupled_periodic_family

SCHEMA="paper03-joint-circle-v3"
FILES=["radius.csv","curvature.csv","shear.csv","bifurcation.csv","twist.csv","invariants.csv","scalars.json","twist_kinetic.csv","shape_minima.csv","cross_kinetic.csv","coupled_action.csv"]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def recompute(output):
    out=Path(output)
    out.mkdir(parents=True,exist_ok=False)
    def csvfile(name,columns,rows):
        with (out/name).open("w",encoding="utf-8",newline="") as f:
            writer=csv.writer(f);writer.writerow(columns);writer.writerows(rows)
    rs=np.linspace(.2,5,241)
    csvfile("radius.csv",["r","MX","AX","VT"],
            zip(rs,potential(rs),potential(rs,V=0),potential(rs,eta=0)))
    rows=[]
    for r in rs:
        d=cartan(r=r)
        rows.append((r,d["P"],d["B"]))
    csvfile("curvature.csv",["r","P","B"],rows)
    eps=np.linspace(-.55,1.75,461)
    h=shear_curvatures(eps);hl=shear_curvatures(eps,path="linear_normalized")
    csvfile("shear.csv",["epsilon","H0","H2","H1","H0_linear","H2_linear","H1_linear"],
            zip(eps,*h,*hl))
    epsb=np.linspace(-.55,-.15,401); qb=shear_branch(epsb)
    delta=[potential(3,M=shape_matrix(ep,shear=np.exp(q/np.sqrt(2))-1))
           -potential(3,M=shape_matrix(ep)) for ep,q in zip(epsb,qb)]
    csvfile("bifurcation.csv",["epsilon","q2_positive","delta_U"],zip(epsb,qb,delta))
    rows=[]
    for w in np.linspace(-1,1,201):
        u=float(potential(3,omega=(w,0,0)))
        d=cartan(omega=(w,0,0),delta=(.4,-.2,.3))
        rows.append((w,u,-d["volume"]*d["R_scalar"]/2))
    csvfile("twist.csv",["omega1","U","U_rotated_frame"],rows)
    rows=[]
    T=shear_basis()[0];K=3*np.pi**2
    for t in np.linspace(-.35,.35,281):
        h=t*T;i2=float(np.trace(h@h));i3=float(np.trace(h@h@h))
        exact=float(potential(3,M=exp_symmetric(h))-potential(3))
        quartic=K*(4*i2+40*i3/3+14*i2*i2/3)
        rows.append((t,exact,quartic))
    csvfile("invariants.csv",["q1","delta_U","quartic"],rows)
    rows=[]
    for rho in np.linspace(0,1,201):
        static,derivative,G=twist_bulk((rho,0,0),(0,2*np.pi*rho,0))
        rows.append((rho,static,derivative,G[1,1]))
    csvfile("twist_kinetic.csv",["rho","static","rotating_derivative","G_transverse"],rows)
    epsm=np.linspace(-.7,.6,651);a=(1+epsm)**(1/3)
    rhom,fmin=shape_minimum(epsm)
    csvfile("shape_minima.csv",["epsilon","rho","F_axis","F_min","F_global"],
            zip(epsm,rhom,a**-8-4/a**2,fmin,np.full_like(epsm,-3.)))
    ws=np.linspace(-.8,.8,321);cross=diagonal_cross_coefficients(ws)
    csvfile("cross_kinetic.csv",["w","G_rw","G_sw","G_rs","G_rL","G_sL"],
            zip(ws,*(cross[name] for name in ["G_rw","G_sw","G_rs","G_rL","G_sL"])))
    coupled=coupled_periodic_family()
    columns=["n","static","derivative","total","delta_static","delta_total"]
    csvfile("coupled_action.csv",columns,zip(*(coupled[name] for name in columns)))
    scalars={}
    for name,eta,V in [("MX",1,1),("AX",1,0),("VT",0,1)]:
        J,H=scalar_data(eta=eta,V=V)
        scalars[name]={"coordinates":["r","L","eta","V"],"point":[3,1,eta,V],
                       "gradient":J.tolist(),"hessian":H.tolist(),
                       "hessian_eigenvalues":np.linalg.eigvalsh(H).tolist()}
    (out/"scalars.json").write_text(json.dumps(scalars,indent=2)+"\n",encoding="utf-8")
    manifest={"schema":SCHEMA,"version":"v3","parameters":{"r":3,"L":1,"kappa":1,"eta":1,"V":1,
                 "theta":"arbitrary constant; complete NY integrates to zero for periodic profiles","alpha":0},
              "geometry":{"circle_period":1,"round_radius":"r/2","volume":"pi**2*L*r**3/4",
                 "torsion":"geometric torsion T as in the manuscript",
                 "shape_path":"M0 exp(q T_A) unless explicitly labelled linear_normalized",
                 "shape_background":"M0=diag(a,a,a^-2); a=(1+epsilon)^(1/3)"},
              "grids":{"radius":[.2,5,241],"shear":[-.55,1.75,461],
                       "bifurcation":[-.55,-.15,401],"twist":[-1,1,201],"invariants":[-.35,.35,281],
                       "twist_kinetic":[0,1,201],"shape_minima":[-.7,.6,651],
                       "cross_kinetic":[-.8,.8,321],"coupled_frequency":[1,16,16],
                       "coupled_periodic_samples":512},
              "circle_variations":{"twist":"fixed r,L,eta,V; omega=rho(cos(2*pi*y),sin(2*pi*y),0)",
                  "scalar":"round, omega=0; periodic r,L,eta,V; bulk after total derivatives",
                  "joint":"r,L,eta,V,M,omega periodic; det M=1; equation M6",
                  "coupled_example":{"amplitude":.5,"shape":"diag(sqrt(z),1/sqrt(z),1); z=1-L^2*w^2/(2*r^2)",
                      "twist":"(A*cos(2*pi*n*y),0,0)","C_A":coupled["C_A"],"round_static":coupled["round_static"]},
                  "interpretation":"Static potential, restricted circle actions and joint Euclidean bulk action"},
              "model_sha256":sha(Path(__file__).with_name("model.py")),
              "circle_sha256":sha(Path(__file__).with_name("circle.py")),
              "data_sha256":{name:sha(out/name) for name in FILES}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return out

def load_data(source):
    src=Path(source);manifest=json.loads((src/"manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema")!=SCHEMA or set(manifest.get("data_sha256",{}))!=set(FILES):
        raise ValueError("Unsupported dataset schema.")
    for name,expected in manifest["data_sha256"].items():
        if sha(src/name)!=expected:raise ValueError(f"Dataset hash mismatch: {name}")
    return {name:np.genfromtxt(src/name,delimiter=",",names=True) for name in FILES if name.endswith(".csv")},json.loads((src/"scalars.json").read_text(encoding="utf-8"))
