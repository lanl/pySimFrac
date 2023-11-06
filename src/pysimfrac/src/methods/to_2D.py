import numpy as np



def generate_2D(self, ind=None, mean_aperture=28):
    
    """Generate a 2D representation of the fracture.

    Parameters
    -------------
        ind : int, optional
            Index along the x-axis to select the slice for 2D representation. If not specified, defaults to the middle of the array.
        mean_aperture : float, optional
            Target mean aperture value to set for the 2D fracture representation. Default is 28.

    Notes
    ---------
        This method generates a 2D representation of the fracture by tiling a selected slice across the x-axis.
        The mean aperture of the fracture is then adjusted to the specified value using the `set_mean_aperture` method.
        The selected slice for the 2D representation can be specified by the index; if not, the middle slice is used by default.
        The attributes 'top' and 'bottom' of the object are modified in-place to represent the 2D fracture.
    """
    
    
    if not ind:
        ind = self.nx//2
    
    self.top = np.tile( self.top[ind], (self.top.shape[1], 1))
    self.bottom = np.tile( self.bottom[ind], (self.top.shape[1], 1))
    self.set_mean_aperture(mean_aperture) 




