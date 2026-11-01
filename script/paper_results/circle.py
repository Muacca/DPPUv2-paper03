"""Joint Euclidean circle action for the coframe specified in the paper.

All profiles are spatially homogeneous. The returned bulk density has periodic
total derivatives removed; it is not the local curvature density. kappa is
positive, theta is constant, and the curvature-squared coefficient is zero.
"""
import numpy as np
from .model import potential


def _array(value, shape, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real.")
    a=np.asarray(value,dtype=float)
    if a.shape!=shape or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be a finite array of shape {shape}.")
    return a


def joint_bulk(r=3.,L=1.,M=None,omega=(0.,0.,0.),*,r_prime=0.,L_prime=0.,
               M_prime=None,omega_prime=(0.,0.,0.),eta=1.,V=1.,kappa=1.):
    """Evaluate equations (M4)-(M6) at one circle coordinate.

    M has determinant one and M_prime satisfies tr(M^-1 M_prime)=0.
    eta and V may be the instantaneous values of periodic scalar profiles;
    their derivatives and coframe rotations contribute only boundary currents.
    The velocity-linear contribution includes the non-Killing spatial shift.
    A positive value on one restricted family is not a global action bound.
    """
    scalars=np.asarray([r,L,r_prime,L_prime,eta,V,kappa])
    if np.iscomplexobj(scalars) or not np.all(np.isfinite(scalars)):
        raise ValueError("Scalar arguments must be finite and real.")
    if r<=0 or L<=0 or kappa<=0:
        raise ValueError("r, L and kappa must be positive.")
    M=_array(np.eye(3) if M is None else M,(3,3),"M")
    Mp=_array(np.zeros((3,3)) if M_prime is None else M_prime,(3,3),"M_prime")
    w=_array(omega,(3,),"omega");wp=_array(omega_prime,(3,),"omega_prime")
    if not np.isclose(np.linalg.det(M),1.,rtol=1e-10,atol=1e-10):
        raise ValueError("M must have determinant one.")
    if not np.isclose(np.trace(np.linalg.solve(M,Mp)),0.,atol=1e-10):
        raise ValueError("M_prime must be tangent to det M=1.")
    B=M.T@M;Bp=Mp.T@M+M.T@Mp
    gamma=r*r*B+L*L*np.outer(w,w)
    gp=2*r*r_prime*B+r*r*Bp+2*L*L_prime*np.outer(w,w)+L*L*(np.outer(wp,w)+np.outer(w,wp))
    inv=np.linalg.inv(gamma);J=L*r**3;N2=J*J/np.linalg.det(gamma)
    beta=L*L*inv@w
    Omega=4*np.array([[0.,-beta[2],beta[1]],[beta[2],0.,-beta[0]],[-beta[1],beta[0],0.]])
    lie=gamma@Omega-Omega@gamma
    X=inv@gp;Y=inv@lie
    coefficient=np.pi**2*J/(32*kappa*kappa*N2)
    quadratic=float(coefficient*(np.trace(X@X)-np.trace(X)**2))
    linear=float(-2*coefficient*np.trace(X@Y))
    static=float(potential(r,L,kappa,eta,V,M,w))
    return {"static":static,"velocity_quadratic":quadratic,"velocity_linear":linear,
            "bulk":static+quadratic+linear,"gamma":gamma,"lapse":float(np.sqrt(N2)),"shift":beta}


def diagonal_cross_coefficients(w,r=3.,L=1.,s=0.,kappa=1.):
    """Velocity cross coefficients (X2), M=diag(exp(s),exp(-s),1)."""
    w=np.asarray(w,dtype=float)
    if r<=0 or L<=0 or kappa<=0 or not np.all(np.isfinite(w)) or not np.isfinite(s):
        raise ValueError("Use positive scales and finite coordinates.")
    c=np.pi**2/(4*kappa*kappa*np.exp(2*s))
    return {"G_rw":-2*c*L*w,"G_sw":c*L*r*w,"G_rs":c*L*w*w,
            "G_rL":-2*c*w*w,"G_sL":c*r*w*w}


def coupled_periodic_family(frequencies=range(1,17),amplitude=.5,*,r=3.,L=1.,
                            kappa=1.,eta=1.,V=1.,samples=512):
    """Fixed-grid periodic quadrature of the exact family (X4).

    Integrate one repeated period and use the exact n^2 derivative scaling.
    No differentiation of tabulated data or high-frequency grid extrapolation
    is used. The returned coefficient C_A is positive for nonzero amplitude.
    """
    ns=np.asarray(list(frequencies),dtype=float)
    if ns.ndim!=1 or len(ns)==0 or not np.all(np.isfinite(ns)) or np.any(ns<1) or np.any(ns!=np.floor(ns)):
        raise ValueError("frequencies must contain positive integers.")
    if r<=0 or L<=0 or kappa<=0 or not np.isfinite(amplitude) or not 0<amplitude<np.sqrt(2)*r/L:
        raise ValueError("Use positive scales and 0 < amplitude < sqrt(2)*r/L.")
    if not isinstance(samples,(int,np.integer)) or samples<64:
        raise ValueError("samples must be an integer at least 64.")
    phase=2*np.pi*np.arange(samples)/samples
    w=amplitude*np.cos(phase);wp=-2*np.pi*amplitude*np.sin(phase)
    z=1-L*L*w*w/(2*r*r);t=np.sqrt(z)
    static=float(np.mean([potential(r,L,kappa,eta,V,np.diag([a,1/a,1.]),(b,0,0)) for a,b in zip(t,w)]))
    C=float(np.mean(np.pi**2*L**3*w*w*wp*wp/(16*kappa*kappa*r*z*z)))
    baseline=float(potential(r,L,kappa,eta,V))
    derivative=-C*ns*ns;total=static+derivative
    return {"n":ns.astype(int),"static":np.full_like(ns,static),"derivative":derivative,
            "total":total,"delta_static":np.full_like(ns,static-baseline),"delta_total":total-baseline,
            "C_A":C,"round_static":baseline}
