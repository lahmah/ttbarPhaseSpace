import numpy as np
from copy import copy

class Vec4D(object):
    def __init__(self, arr=None):
        if arr is None:
            self._arr = np.zeros(4)
        else:
            arr = np.asarray(arr)
           # if arr.shape != (4,):
           #     raise TypeError('Wrong array size! Must have 4 entries.')

            self._arr = arr

    def __add__(self, rhs):
        if isinstance(rhs, self.__class__):
            return self.__class__(self._arr + rhs._arr)
        elif isinstance(rhs, int) or isinstance(rhs, float):
            return self.__class__(self._arr + rhs)
        else:
            return NotImplemented

    def __sub__(self, rhs):
        if isinstance(rhs, self.__class__):
            return self.__class__(self._arr - rhs._arr)
        elif isinstance(rhs, int) or isinstance(rhs, float):
            return self.__class__(self._arr - rhs)
        else:
            return NotImplemented

    def __mul__(self, rhs):
        if isinstance(rhs, self.__class__):
            mk = np.array([-1,1,1,1])
            if len(self._arr.shape) == 2:
                mk = mk[:,None]
            return np.sum(-self._arr*rhs._arr*mk,0)
        elif isinstance(rhs, int) or isinstance(rhs, float):
            return self.__class__(self._arr * rhs)
        else:
            return NotImplemented

    def __rmul__(self, lhs):
        if isinstance(lhs, int) or isinstance(lhs, float):
            return self.__class__(lhs * self._arr)
        else:
            return NotImplemented

    def __getitem__(self, key):
        if hasattr(key,"__len__"):
            if key[0] == slice(None,None,None) and isinstance(key[1],int):
                return Mom4D(self._arr[key])
        
        return self._arr[key]

    def __setitem__(self, key, value):
        self._arr[key] = value

    @property
    def vec3d(self):
        return self[1:]

    @vec3d.setter
    def vec3d(self, arr):
        arr = np.asarray(arr)
        if arr.shape != (3,):
            raise TypeError('Wrong array size! Must have 3 entries.')

        self[1:] = arr

    def __str__(self):
        return str(self._arr)

    def __repr__(self):
        return str(self._arr)
class Mom4D(Vec4D):
    @property
    def mom3d(self):
        return self.vec3d

    @mom3d.setter
    def mom3d(self, p):
        self.vec3d = p

    @property
    def E(self):
        return self[0]

    @E.setter
    def E(self, E):
        self[0] = E

    @property
    def m(self):
        m2 = self*self
        if isinstance(m2,float):
            m2 = np.array([m2])
        m2[np.logical_and(np.isclose(0,m2,atol=1.e-7),m2<0)]=0
        if (m2 < 0).any():
               raise ValueError('Negative Mass!')

        if len(m2)==1:
            return np.sqrt(m2[0])
        return np.sqrt(m2)

    @property
    def pT(self):
        return np.sqrt(self[1]*self[1]+self[2]*self[2])

# angle between two 4-momenta
def angle(p, q):
    if not (isinstance(p, Mom4D) and isinstance(q, Mom4D)):
        raise TypeError("Arguments need to be of type 'Mom4D'")

    cos_angle = (p.mom3d.dot(q.mom3d))/(p.E*q.E)  
    
    return cos_angle
