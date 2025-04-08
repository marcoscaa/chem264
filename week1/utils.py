import numpy as np

def make_pos(N, R):
    """
    Generates an equally spaced array with N numbers, symmetric around 0,
    and with a spacing of R.
   
    Args:
      N: The number of H atoms.
      R: The internuclear spacing.
   
    Returns:
      A numpy array of equally spaced numbers symmetric around 0.
    """
   
    if N % 2 == 0:  # Even number of elements
      half_N = N // 2
      array = np.arange(-R * (half_N - 0.5), R * half_N, R)
    else:  # Odd number of elements
      half_N = N // 2
      array = np.arange(-R * half_N, R * (half_N + 1), R)
   
    return array

def discretize_space(pos,dx=0.01):
    xmin=np.min(pos)-5.0
    xmax=np.max(pos)+5.0
    return np.arange(xmin,xmax,dx)
