from typing import List, Union, Optional
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline, PPoly
from scipy.interpolate import splprep, splev
from scipy import integrate
import math
import enum
import copy
from scipy.optimize import minimize, brentq
from .resample import clean_curve_points
import numpy.typing as npt 
from .order_pts import order_points_by_spline_arc_length

def convert_to_ndarray(t):
    """
        converts a scalar or list to numpy array 
    """

    if type(t) is not np.ndarray and type(t) is not list: # Scalar
        t = np.array([t],dtype=float)
    elif (type(t) is list):
        t = np.array(t,dtype=float)
    return t


class spline_type(enum.Enum):
    pchip = 1
    spline = 2 # Cubic spline 

def spline_intersect(spline1,spline2,xmin,xmax):
    def intersect_test(x0):
        y1 = spline1(x0)
        y2 = spline2(x0)
        return abs(y1-y2)
    bnds = ((xmin,xmax),)
    x0 = (xmin+xmax)/2.0
    res = minimize(intersect_test,x0=(x0),bounds=bnds,tol=1E-5)
    if (res.success):
        return res.x[0]
    else:
        return -1

def pspline_intersect(pspline1,pspline2,tmin,tmax):
    def intersect_test(t0):
        pt1,_ = pspline1(t0[0])
        pt2,_ = pspline2(t0[1])
        dx = pt2[0,0]-pt1[0,0]
        dy = pt2[0,1]-pt1[0,1]
        return math.sqrt(dx*dx + dy*dy)

    bnds = ((tmin,tmax),(tmin,tmax),)
    t0 = (tmin+tmax)/2.0
    res = minimize(intersect_test,x0=(t0,t0),bounds=bnds,tol=1E-4)
    if (res.success and res.fun<0.005):
        return res.x
    else:
        return res.x*0-1 # return -1 for no intersect


 
class pspline:
    """fits a spline on a curve and allows you to extract points with a percentage as input"""
    pxy:npt.NDArray
    ndim:int
    chordlen:float
    cumarclen:npt.NDArray
    spl:List[Union[PchipInterpolator,CubicSpline]]
    spld:List[float]
    polyarray:npt.NDArray
    totalsplinelength:float
    seglen:npt.NDArray
    
    
    def __segkernel__(self,y,t):
        t = convert_to_ndarray(t)
        val = np.zeros(len(t))
        for k in range(0,self.ndim):
            val += np.polyval(self.polyarray[k,:],t)**2
        val = np.sqrt(val)
        return val
    
    def __init__(self,px:Union[npt.NDArray,List[float]],py:Union[npt.NDArray,List[float]],pz:Optional[Union[npt.NDArray,List[float]]]=[],method=spline_type.pchip):
        """Creates a percentage spline

        Args:
            px (Union[npt.NDArray,List[float]]): list of x values
            py (Union[npt.NDArray,List[float]]): list of y values 
            pz (Union[npt.NDArray,List[float]], optional): Optional, list of z vlaues. Defaults to [].
            method (spline_type, optional): Type of spline to define. Defaults to spline_type.pchip.
        """
        px = convert_to_ndarray(px)
        py = convert_to_ndarray(py)
                
        coordinates = (px, py) if pz is None or len(pz) == 0 else (px, py, pz)
        self.pxy = clean_curve_points(np.stack(coordinates, axis=1))
        self.ndim = self.pxy.shape[1]
        self.chordlen = np.linalg.norm(np.diff(self.pxy, axis=0), axis=1)
        self.n = self.pxy.shape[0]
        self.chordlen = self.chordlen/np.sum(self.chordlen)
        self.cumarclen = np.append(np.zeros(1),np.cumsum(self.chordlen))

        # compute parametric splines
        self.spl = []
        self.spld = []
        diffarray = np.array([[3,0,0],[0,2,0],[0,0,1],[0,0,0]])
        for i in range(0,self.ndim):
            if (method==spline_type.pchip):
                self.spl.append(PchipInterpolator(self.cumarclen,self.pxy[:,i]))
            if (method==spline_type.spline):
                self.spl.append(CubicSpline(self.cumarclen,self.pxy[:,i]))
                # nc = numel(self.spl[i].coefs)
                # if nc < 4:
                #     # just pretend it has cubic segments
                #     self.spl[i].coefs = np.stack(np.zeros(1,4-nc),self.spl[i].coefs)
                #     self.spl[i].order = 4
           
            # and now differentiate them
            xp = copy.deepcopy(self.spl[i])
            xp.c = np.transpose(np.matmul(np.transpose(xp.c),diffarray))
            # self.xp.order = 3
            self.spld.append(xp)

        '''
            Generate the total arclength along the curve
            by integrating each segment and summing the
            results. The integration scheme does its job
            using an ode solver.
        '''
        
        # polyarray here contains the derivative polynomials
        # for each spline in a given segment
        self.polyarray = np.zeros((self.ndim,3))
        self.seglen = np.zeros((self.n-1,1))
        _,pieces = self.spl[0].c.shape
        for i in range(pieces):
            # extract polynomials for the derivatives
            for j in range(self.ndim):
                self.polyarray[j,:] = self.spld[j].c[:,i] # Grab the coefficients
            
            '''
                Integrate the arclength for the i'th segment
                using ode45 for the integral. I could have
                done self part with quad too, but then it
                would not have been perfectly (numerically)
                consistent with the next operation in self tool.
            '''
            self.seglen[i] = integrate.quad(lambda u: float(self.__segkernel__(0, u)[0]), 0, self.chordlen[i], epsabs=1e-10)[0]
            # r = integrate.ode(self.__segkernel__).set_integrator('vode',method='bdf') 
            # r.set_initial_value(0,0)    
            
            # self.seglen[i] = r.integrate(1)
        
        # and normalize the segments to have unit total length
        self.totalsplinelength = np.sum(self.seglen)
        self.cumseglen = np.append(np.zeros(1),np.cumsum(self.seglen))

    def __call__(self, t):
        """:no-index:"""
        return self.get_point(t)

    def get_curve_len(self,t):
        return self.totalsplinelength*t
    
    def get_point(self,t:Union[float,List[float],npt.NDArray]) -> npt.NDArray:
        """Gets the points anywhere along the spline 

        Args:
            t (Union[float,List[float]]): array or float from 0 to 1

        Returns:
            Tuple containing:
                
                (npt.NDArray): x,y,z value of the spline 
                (npt.NDArray): dudt - derivative with time
        """
    
        fractions = np.atleast_1d(np.asarray(t, dtype=float))
        if fractions.ndim != 1 or not np.isfinite(fractions).all() or np.any((fractions < 0) | (fractions > 1)):
            raise ValueError("Arc fractions must be finite and within [0, 1]")
        parameters = np.empty_like(fractions)
        # Invert integrated spline speed within each segment. Preserve query order.
        for i, fraction in enumerate(fractions):
            if fraction == 0 or fraction == 1:
                parameters[i] = fraction
                continue
            target = fraction * self.totalsplinelength
            segment = int(np.clip(np.searchsorted(self.cumseglen, target, side="right") - 1, 0, self.n - 2))
            lo, hi = self.cumarclen[segment:segment + 2]
            remaining = target - self.cumseglen[segment]
            def speed(u):
                return float(np.linalg.norm([spline(u) for spline in self.spld]))
            length = integrate.quad(speed, lo, hi, epsabs=1e-10)[0]
            remaining = np.clip(remaining, 0., length)
            parameters[i] = brentq(lambda u: integrate.quad(speed, lo, u, epsabs=1e-10)[0] - remaining, lo, hi)
        points = np.column_stack([spline(parameters) for spline in self.spl])
        derivatives = np.column_stack([spline(parameters) for spline in self.spld])
        return points, derivatives
