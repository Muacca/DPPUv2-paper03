"""Euclidean EC + complete-NY geometry and specified circle variations.

Parameters eta and V are the manuscript's geometric torsion coefficients.
The circle has period one and circumference L; the round
SU(2) radius is r/2. Circle-derivative functions concern the stated restricted
Euclidean families defined in the manuscript.
"""
import itertools
import numpy as np

def shear_basis():
    return np.array([
        np.diag([1.,1.,-2.])/np.sqrt(6),
        np.diag([1.,-1.,0.])/np.sqrt(2),
        [[0,1/np.sqrt(2),0],[1/np.sqrt(2),0,0],[0,0,0]],
        [[0,0,1/np.sqrt(2)],[0,0,0],[1/np.sqrt(2),0,0]],
        [[0,0,0],[0,0,1/np.sqrt(2)],[0,1/np.sqrt(2),0]],
    ])

def exp_symmetric(h):
    values,vectors=np.linalg.eigh(np.asarray(h,dtype=float))
    return (vectors*np.exp(values))@vectors.T

def shape_matrix(epsilon=0., shear=0., q=(0.,0.,0.)):
    if epsilon <= -1 or shear <= -1:
        raise ValueError("epsilon and shear must exceed -1.")
    a=(1+epsilon)**(1/3)
    M=np.diag([a*(1+shear),a/(1+shear),a**-2])
    for x,T in zip(q,shear_basis()[2:]):
        M=M@exp_symmetric(x*T)
    return M

def mixing_rotation(delta=(0.,0.,0.)):
    O=np.eye(4)
    for i,d in enumerate(delta):
        c=1/np.sqrt(1+d*d); s=d*c
        R=np.eye(4);R[i,i]=R[3,3]=c;R[i,3]=s;R[3,i]=-s
        O=O@R
    return O

def epsilon_tensor(n):
    E=np.zeros((n,)*n)
    for p in itertools.permutations(range(n)):
        E[p]=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
    return E

EPS3=epsilon_tensor(3)
EPS4=epsilon_tensor(4)

def cartan(r=3., L=1., eta=1., V=1., M=None, omega=(0.,0.,0.), delta=(0.,0.,0.)):
    """Construct coframe, connection, curvature and densities directly.

    SO(4) mixing rotates the coframe. The torsion ansatz is specified in the
    resulting frame. Keeping its components fixed is not a simultaneous gauge
    transformation of all fields; the scalar action nevertheless has no
    dependence on this orientation on the stated unimodular family.
    """
    if r<=0 or L<=0: raise ValueError("r and L must be positive.")
    M=np.eye(3) if M is None else np.asarray(M,dtype=float)
    if M.shape!=(3,3) or not np.isclose(np.linalg.det(M),1.,rtol=1e-10,atol=1e-10):
        raise ValueError("M must be a real 3x3 matrix with determinant one.")
    A=np.zeros((4,4));A[:3,:3]=r*M;A[3,:3]=L*np.asarray(omega);A[3,3]=L
    A=mixing_rotation(delta)@A
    inv=np.linalg.inv(A)
    base=np.zeros((4,4,4));base[:3,:3,:3]=4*EPS3
    C=np.einsum("ad,eb,fc,def->abc",A,inv,inv,base)
    GammaLC=.5*(C+C.transpose(2,1,0)-C.transpose(1,0,2))
    T=np.zeros((4,4,4));T[:3,:3,:3]=2*eta/r*EPS3
    for a,b,c in itertools.product(range(4),repeat=3):
        T[a,b,c]+=V/3*((a==c)*(b==3)-(a==b)*(c==3))
    K=-.5*(T+T.transpose(2,0,1)-T.transpose(1,2,0))
    Gamma=GammaLC+K
    def curvature(G):
        return (np.einsum("aec,ebd->abcd",G,G)-np.einsum("aed,ebc->abcd",G,G)
                +np.einsum("abe,ecd->abcd",G,C))
    R=curvature(Gamma); RLC=curvature(GammaLC)
    TT=.25*np.einsum("abcd,eab,ecd->",EPS4,T,T)
    REE=.5*np.einsum("abcd,abcd->",EPS4,R)
    P=.25*np.einsum("cdef,abcd,abef->",EPS4,R,R)
    B=sum(R[0,1,c,d]*R[2,3,c,d]-R[0,2,c,d]*R[1,3,c,d]+R[1,2,c,d]*R[0,3,c,d]
          for c in range(4) for d in range(c+1,4))
    T_reconstructed=C+Gamma.transpose(0,2,1)-Gamma
    return {"coframe":A,"C":C,"Gamma":Gamma,"T":T,"R":R,"R_LC":RLC,
            "R_scalar":float(np.einsum("abab->",R)),
            "R_scalar_LC":float(np.einsum("abab->",RLC)),
            "TT":float(TT),"REE":float(REE),"NY":float(TT-REE),
            "P":float(P),"B":float(B),"torsion_residual":float(np.max(np.abs(T_reconstructed-T))),
            "volume":float(np.pi**2*L*r**3/4)}

def potential(r, L=1., kappa=1., eta=1., V=1., M=None, omega=(0.,0.,0.)):
    """Closed form of minus the integrated static EC + complete-NY action."""
    r=np.asarray(r,dtype=float)
    if np.any(r<=0) or L<=0 or kappa<=0: raise ValueError("r,L,kappa must be positive.")
    M=np.eye(3) if M is None else np.asarray(M,dtype=float)
    if not np.isclose(np.linalg.det(M),1.,rtol=1e-10,atol=1e-10):
        raise ValueError("M must have determinant one.")
    S=M@M.T; w=np.asarray(omega,dtype=float)
    return np.pi**2*L*r/kappa**2*(2*np.trace(S@S)-np.trace(S)**2
            +3*eta**2/4+V**2*r**2/12+L**2/r**2*np.dot(M@w,M@w))

def shear_curvatures(epsilon,r=3.,L=1.,kappa=1.,path="exponential"):
    """Second derivatives at M0=diag(a,a,a^-2), r and torsion fixed.

    Exponential path: M=M0 exp(q T_A).
    Linear-normalized path: M=M0 (I+q T_A)/det(I+q T_A)^(1/3).
    The two agree at the stationary round shape, not at a generic shape.
    """
    ep=np.asarray(epsilon,dtype=float)
    if np.any(ep<=-1):raise ValueError("epsilon must exceed -1.")
    a=(1+ep)**(1/3); K=np.pi**2*L*r/kappa**2
    if path=="exponential":
        return np.array([8*K*(4-a**6)/(3*a**8),
                         8*K*(2*a**6-1)/a**2,4*K*(a**6+1)/a**8])
    if path=="linear_normalized":
        return np.array([4*K*(7-a**6)/(3*a**8),
                         4*K*(4*a**6-1)*(3*a**6-1)/(3*a**8),
                         2*K*(7*a**6+5)/(3*a**8)])
    raise ValueError("Unknown path.")

def shear_branch(epsilon):
    """Positive q2 minimum at fixed epsilon; negative branch is its opposite."""
    a6=(1+np.asarray(epsilon,dtype=float))**2
    return np.arccosh(np.maximum(1.,1/(2*a6)))/np.sqrt(2)

def scalar_data(r=3.,L=1.,eta=1.,V=1.,kappa=1.):
    """Gradient and Hessian in geometric coordinates (r,L,eta,V), round, omega=0."""
    p=np.pi**2/kappa**2
    J=p*np.array([L*(3*V*V*r*r+9*eta*eta-36)/12,
                  r*(V*V*r*r+9*eta*eta-36)/12,1.5*L*eta*r,L*V*r**3/6])
    H=p*np.array([[L*V*V*r/2,(3*V*V*r*r+9*eta*eta-36)/12,1.5*L*eta,L*V*r*r/2],
                  [(3*V*V*r*r+9*eta*eta-36)/12,0,1.5*eta*r,V*r**3/6],
                  [1.5*L*eta,1.5*eta*r,1.5*L*r,0],
                  [L*V*r*r/2,V*r**3/6,0,L*r**3/6]])
    return J,H

def shape_minimum(epsilon):
    """Minimize the full (q2,q3) doublet at fixed axial control epsilon.

    Return amplitude rho and dimensionless shape minimum F. All-volume-
    preserving shape minimization additionally varies epsilon and gives F=-3.
    """
    ep=np.asarray(epsilon,dtype=float)
    if np.any(ep<=-1):raise ValueError("epsilon must exceed -1.")
    a=(1+ep)**(1/3)
    rho=shear_branch(ep)
    value=np.where(a**6<.5,-4*a**4,a**-8-4*a**-2)
    return rho,value

def twist_bulk(omega,omega_prime,r=3.,L=1.,kappa=1.):
    """Exact circle-action density above zero twist, fixed scalar controls.

    The round spatial coframe has arbitrary real vector omega(y); the complete
    action is integrated over a period-one circle. Periodic total derivatives
    have been removed. The quartic term is in amplitude, not derivative order.
    Return static density, derivative density and the velocity Hessian.
    """
    if min(r,L,kappa)<=0:raise ValueError("r,L,kappa must be positive.")
    w=np.asarray(omega,dtype=float);p=np.asarray(omega_prime,dtype=float)
    if w.shape!=(3,) or p.shape!=(3,) or not np.all(np.isfinite([w,p])):
        raise ValueError("omega and omega_prime must be finite three-vectors.")
    c=np.pi**2*L**3/(kappa**2*r);w2=np.dot(w,w)
    Q=np.dot(np.cross(w,p),np.cross(w,p))
    return c*w2,c*Q/16,c*(w2*np.eye(3)-np.outer(w,w))/8

def scalar_bulk(r,r_prime,L=1.,eta=1.,V=1.,kappa=1.):
    """Round scalar circle-action density after periodic total derivatives.

    r,L,eta,V may depend on the circle coordinate. Only r_prime remains in
    the bulk derivative term; L here is the local circle lapse.
    """
    if min(r,L,kappa)<=0:raise ValueError("r,L,kappa must be positive.")
    return float(potential(r,L=L,kappa=kappa,eta=eta,V=V)
                 -3*np.pi**2*r*r_prime**2/(4*kappa**2*L))
