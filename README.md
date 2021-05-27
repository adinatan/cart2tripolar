[![View cart2tripolar on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/65217-cart2tripolar)
# cart2tripolar

 Cartesian to polar pixel transform - with compact representation
  The code starts with a Cartesian image of square pixel size dx*dy=1, and
  outputs a polar image of polar pixels that has a similar amount of
  information in polar coordinates dr*dtheta~=1.  The signal in each polar
  pixel is determined via its fractional overlap with the four surrounding
  Cartesian pixels. The result is a triangular polar representation,
 because the number of polar pixels per radius per angle increase with
  radius, i.e. the larger the radius the higher the angular resolution.
  The code was originally used in the context of analyzing images with 
 symmetry around a quadrant so that functionality was kept. 
 The code support NaN valued pixels for masking.

See tripol_der.pdf for more details.

  ![Fig1](https://github.com/adinatan/cart2tripolar/blob/master/fig2.png)
