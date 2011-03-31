// pixel : 20mas
// lambda : 1.644
// obs : 0
//

/*
 * pray_core.i
 *
 * This file is part of the pray package
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 * Initial release D. Gratadour, Aug 2009.
 */
require,"yaodh.i";

plug_in,"pray_core";

struct pray_struct
{
  pointer images;
  pointer ftobject;
  pointer variance;
  pointer fft_ws;
  pointer norm;
  pointer nzer;
  pointer alt;
  pointer mircube;
  pointer ipupil;
  pointer def;
  pointer defPup;
  pointer ishifts;
  pointer jshifts;
  pointer xshifts;
  pointer yshifts;
  pointer deltaFoc;
  pointer poffset;
  float   teldiam;
  float   cobs;
  float   threshold;
  float   lambda;
  float   nphe;
  long    scale;
  long    pupd;
  long    size;
  long    _n;
  long    _n1;
  long    _n2;
  long    nbiter;
  long    neval;
  long    dmgrid;
  long    fit_starpos;
  long    diff_tt;
  long    fit_object;
};

/*
 _____                     __   __ _    ___  
|  ___| __ ___  _ __ ___   \ \ / // \  / _ \ 
| |_ | '__/ _ \| '_ ` _ \   \ V // _ \| | | |
|  _|| | | (_) | | | | | |   | |/ ___ \ |_| |
|_|  |_|  \___/|_| |_| |_|   |_/_/   \_\___/ 
                                            

the following routines have been copy-paste from yao
*/

extern _get2dPhase
/* PROTOTYPE
   int _get2dPhase(pointer pscreens, int psnx, int psny, int nscreens, pointer outphase, int phnx, int phny, pointer ishifts, pointer xshifts, pointer jshifts, pointer yshifts)
*/

//----------------------------------------------------

func make_pupil(dim,pupd,xc=,yc=,real=,cobs=)

  /* DOCUMENT func make_pupil(dim,pupd,xc=,yc=,real=)
   */
{
  if (real == 1) {
    pup = exp(-(mydist(dim,xc=xc,yc=yc)/(pupd/2.))^60.)^0.69314;
  } else {
    //    xc;yc;info,xc;
    //    tv,dist(dim,xc=xc,yc=yc);pause,2000;
    pup = mydist(dim,xc=xc,yc=yc) < (pupd+1.)/2.;
  }
  if (cobs!=[]) {
    if (real == 1) {
      pup -= exp(-(mydist(dim,xc=xc,yc=yc)/(pupd*cobs/2.))^60.)^0.69314;
    } else {
      pup -= mydist(dim,xc=xc,yc=yc) < (pupd*cobs+1.)/2.;
    }
  }
    
  return pup;
}

extern _dist
/* PROTOTYPE
   int _dist(pointer dptr, long dimx, long dimy, float xc, float yc)
*/

func mydist(dim,xc=,yc=)
/* DOCUMENT func dist(dim,xc=,yc=)
 * Return an array which elements are the distance to (xc,yc). xc and
 * yc can be omitted, in which case they are defaulted to size/2+1.
 * F.Rigaut, 2003/12/10.
 * SEE ALSO:
*/

{
  dim = long(dim);

  if (is_scalar(dim)) dim=[2,dim,dim];
  if ((is_vector(dim)) && (dim(1)!=2))
    error,"Dist only deals with 2D square arrays";

  d = array(float,dim);

  if (xc!=[]) {xc = float(xc-1.);} else {xc = float(dim(2)/2);}
  if (yc!=[]) {yc = float(yc-1.);} else {yc = float(dim(3)/2);}

  res = _dist(&d,dim(2),dim(3),xc,yc);

  return d;
}

//---------------------------------------------------------

func zernumero(zn)
/* DOCUMENT zernumero(zn)
 * Returns the radial degree and the azimuthal number of zernike
 * number zn, according to Noll numbering (Noll, JOSA, 1976)
 * SEE ALSO: prepzernike, zernike
 */
{
  j	= 0;
  for (n=0;n<=100;n++)
   {
    for (m=0;m<=n;m++)
     {
       if (((n-m) % 2) == 0)
       {
        j	= j+1;
        if (j == zn) {return [n,m];}
        if (m != 0)
         {
          j	= j+1;
          if (j == zn) {return [n,m];}
         }
       }
     }
   }
}

//---------------------------------------------------------

func gamma(arg)
/* DOCUMENT gamma(arg)
 * Gamma function.
 * SEE ALSO: gamma.i in yorick/i/
 */
{
  if (arg == 0.) {
    return 1.;
  } else {
    return exp(ln_gamma(arg));
  }
}

//---------------------------------------------------------

func factoriel(arg)
/* DOCUMENT factoriel(arg)
 * Return factoriel of the argument    
 * SEE ALSO:
 */
{
  if (arg == 0) {
    return 1.;
  } else {
    res = 1.;
    for (i=1;i<=arg;i++) res = res*i;
    return res;
   }
}

//---------------------------------------------------------

func prepzernike(size,diameter,xc,yc)
/* DOCUMENT prepzernike(size,diameter,xc,yc)
 * Call this function to set up the geometry for subsequent calls
 * to the zernike function.
 * size : size of the 2d array on which future "zernike" will be returned
 * diameter : diameter of the pupil in pixel in the array
 * xc, yc (optional) : Coordinates (in pixels of the center of the pupil)
 * Example:
 * > prepzernike,128,100
 * > pli,zernike(6)
 * SEE ALSO: zernike,zernike_ext,zernumero
 */
{
  extern zdim,zr,zteta,zmask,zrmod,zmaskmod;

  if (xc == []) {xc = size/2+1;}
  if (yc == []) {yc = size/2+1;}

  radius= (diameter+1.)/2.;
  zdim	= size;
  zr	= mydist(zdim,xc=xc,yc=yc)/radius;
  zmask	= (zr <= 1.);
  zmaskmod = (zr <= 1.2);
  zrmod	= zr*zmaskmod;
  zr	= zr*zmask;
  x	= float(span(1,zdim,zdim)(,-:1:zdim));
  y	= transpose(x);
  zteta	= atan(y-yc,x-xc);
}

//---------------------------------------------------------

func zernike(zn)
/* DOCUMENT zernike(zn)
 * Returns the zernike number zn, defined on a 2D array as per
 * the prepzernike function.
 * These zernikes follow the Noll (JOSA, 1976) numbering and
 * definition (rms of 1 over the pupil)
 * Example:
 * > prepzernike,128,100
 * > pli,zernike(6)
 * SEE ALSO: prepzernike, zernumero
 */
{
  extern zdim,zr,zteta,zmask,zrmod,zmaskmod;

  z	= array(float,zdim,zdim);
  znm	= zernumero(zn) ; n=znm(1) ; m=znm(2);

  for (i=0;i<=(n-m)/2;i++) {
    z = z + (-1.)^i*zr^(n-2.*i)*factoriel(n-i)/
      (factoriel(i)*factoriel((n+m)/2-i)*factoriel((n-m)/2-i));
  }
  if ((zn % 2) == 1) {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*sin(m*zteta);
    }
  } else {
    if (m == 0) {
      z = z*sqrt(n+1.);
    } else {
      z = z*sqrt(2*(n+1.))*cos(m*zteta);
    }
  }

  return z*zmask;
}

//----------------------------------------------------

func make_pzt_dm(size,dim,nxact,cobs,pitch,coupling,xflip=,yflip=,pitchMargin=,unitpervolt=)
  /* DOCUMENT function make_pzt_dm2(dm_structure,disp=)
     the influence functions are in microns per volt.
  */
{
  if (xflip == []) xflip = 0;
  if (yflip == []) yflip = 0;
  if (pitchMargin == []) pitchMargin = 0;
  if (unitpervolt == []) unitpervolt = 1;
  
  cent = size/2+0.5;

  // best parameters, as determined by a multi-dimensional fit
  // (see coupling3.i)
  a=[4.49469,7.25509,-32.1948,17.9493];
  p1 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [2.49456,-0.65952,8.78886,-6.23701];
  p2 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [1.16136,2.97422,-13.2381,20.4395];
  irc = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  /*
    ir  = pitch*1.2;
    ir  = pitch*1.46;
    ir  = pitch*1.65;
    c  = 3.8; p1 = 3.9; p2 = 2.4; ir = pitch*1.20;  // good. no. coupling 8% too low
    c  = 3.8; p1 = 4; p2 = 2.4; ir = pitch*1.65;
    c  = 3.75; p1 = 4.2; p2 = 2.5; ir = pitch*1.25; //ok, coupling=13%
    c  = 3.74; p1 = 3.805; p2 = 2.451; ir = pitch*1.4; //good, coupling=17%
    c  = 4; p1 = 3.84; p2 = 2.5; ir = pitch*1.5; //good, coupling=20%
    c  = 3.74; p1 = 3.805; p2 = 2.451; ir = pitch*1.4; //good, coupling=17%
  */
  ir = irc*pitch;

  bord  = 0;
  cub   = array(float,nxact+bord*2,nxact+bord*2,4);

  // make X and Y indices array:
  xy    = indices(nxact+bord*2);

  // express "centered" coordinate of actuator in pixels:
  xy    = (xy-1.-bord-(nxact-1.)/2.)*pitch;

  // fill cub (X coord  and Y coord):
  cub(,,1) = xy(,,1); cub(,,2) = xy(,,2);

  if (xflip) cub(,,1) = cub(::-1,,1);
  if (yflip) cub(,,2) = cub(,::-1,1);

  dis      = sqrt(cub(,,1)^2.+cub(,,2)^2.);
  
  if (pitchMargin == 0)  pitchMargin = 1.44;

  rad      = ((nxact-1.)/2.+pitchMargin)*pitch; //+1.44 is the margin
  inbigcirc= where(dis < rad);
  // 1 if valid actuator, 0 if not:
  // selection is done after interaction matrix is done
  cub(,,3) = 1;

  // 1 if valid guard ring actuator, 0 if not:
  //cub(,,4) = (dis >= (pupr+extent*pitch)) & (dis < (pupr+(1.+extent)*pitch));
 // I don't use extrapolation actuator anymore.
  cub(,,4) = 0.;

  // converting to array coordinates:
  cub(,,1) = cub(,,1)+cent;
  cub(,,2) = cub(,,2)+cent;

  cub      = cub(*,);
  // cub now has two indices: first one is actuator number (valid or extrap)
  // second one is: 1:Xcoord, 2:Ycoord, 3:valid?, 4:extrapolation actuator?

  // filtering actuators outside of a disk radius = rad (see above)
  cub      = cub(inbigcirc,);

  cubval   = cub(where(cub(,3)),);

  nvalid   = int(sum(cubval(,3)));

  xy    = indices(size);
//  x     = xy(,,2); y = xy(,,1);
  x     = xy(,,1); y = xy(,,2);
  def = array(float,size,size,nvalid);

  extent = pitch*(nxact+2.); // + 1.5 pitch each side      
  _n1 = long(clip(floor(cent-extent/2.),1,));
  _n2 = long(clip(ceil(cent+extent/2.),,dim));
  //x = x(_n1:_n2,_n1:_n2);
  //y = y(_n1:_n2,_n1:_n2);

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  for (i=1;i<=nvalid;i++) {
    tmpx       = clip(abs((x-cubval(i,1))/ir),1e-8,2.);
    tmpy       = clip(abs((y-cubval(i,2))/ir),1e-8,2.);
    tmp        = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*(1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
    def(,,i)   = tmp*(tmpx <= 1.)*(tmpy <= 1.);
  }

  tmp=pitch/abs(ir);
  coupling = 1.- tmp^p1 + c*log(tmp)*tmp^p2;

  fact = unitpervolt/max(def);

  def = float(def*fact);

  return def;
}


func filter_quad(size,dim,centx,centy,basis,foconly=)
{
  if (foconly == []) foconly = 0;
  prepzernike,size,dim,centx,centy;
  valids = where(zernike(1));
  if (foconly) quads = [zernike(2),zernike(3),zernike(4)];
  else quads = [zernike(2),zernike(3),zernike(4),zernike(5),zernike(6)];
  
  matPass = (LUsolve(quads(*,)(valids,)(+,)*quads(*,)(valids,)(+,))(+,)*quads(*,)(valids,)(*,)(,+));
  //matFilt = quads(,,+)*matPas(,+);
  modes = matPass(,+)*basis(*,)(valids,)(+,);

  basis2 = quads(,,+)*modes(+,);

  /*
  matPass = (LUsolve(basis(*,)(valids,)(+,)*basis(*,)(valids,)(+,))(+,)*basis(*,)(valids,)(*,)(,+));
  modes2 = matPass(,+)*quads(*,)(valids,)(+,);

  //res=(dm0(+)*modes2(+,))(+)*modes(+,);

  matFilt = modes2(,+)*modes(+,);

  //res2 = dm0(+)*matFilt(+,);

  matFilt = modes(+,)*modes2(,+);
  basis2 = basis(,,+)*matFilt(,+);

  error;
  */
  return basis-basis2;
}

/*
 ____  _                      ____  _                    _ _         
|  _ \| |__   __ _ ___  ___  |  _ \(_)_   _____ _ __ ___(_) |_ _   _ 
| |_) | '_ \ / _` / __|/ _ \ | | | | \ \ / / _ \ '__/ __| | __| | | |
|  __/| | | | (_| \__ \  __/ | |_| | |\ V /  __/ |  \__ \ | |_| |_| |
|_|   |_| |_|\__,_|___/\___| |____/|_| \_/ \___|_|  |___/_|\__|\__, |
                                                               |___/ 
*/

func prepmirror(size,patchDiam,nactu,couplingFact,dmpitchX,dmpitchY,x0,y0,defpup,deftt,&gradG,orthonorm=)
/* DOCUMENT 
 * prepmirror(size,patchDiam,nactu,couplingFact,dmpitchX,dmpitchY,x0,y0,defpup,deftt,&gradG,orthonorm=)
 * This routine returns the influence functions in the grid model
 *
 */ 
{
  if (orthonorm == []) orthonorm = 0;

  nactu = long(nactu(1));
  pup = defpup;
  x = span(-1,1,nactu)(,-:1:nactu);
  y = transpose(x);
  k=0;
  
  def = array(0.0, size, size, nactu*nactu);
  gradG = array(0.0, size, size, nactu*nactu,4);
  Ax = -log(couplingFact)/((2./(nactu-1))^2);
  Ay = -log(couplingFact)/((2./(nactu-1))^2);
  
  for(i=1;i<=nactu;i++) {
    for(j=1;j<=nactu;j++) {
      r = sqrt(x(i,j)^2 + y(i,j)^2 );
      if( r<1.1 ) {
        k++;
        xc = x0 + j*dmpitchX;
        yc = y0 + i*dmpitchY;
        def(,,k) = mygauss2(size,xc,yc,Ax/2.82843,Ay/2.82843,1.,0.,0.,grad,deriv=1);
        gradG(,,k,1) = grad(,,1); // grad x0
        gradG(,,k,2) = grad(,,2); // grad y0
        gradG(,,k,3) = grad(,,1)*j; // grad pitchX
        gradG(,,k,4) = grad(,,2)*i; // grad pitchY
      }
    }
  }
  def = pup * def(,,:k);
  gradG = pup * gradG(,,:k,);

  nmodes = dimsof(def)(4);
  norm = sum(deftt(,,2)^2);

  for (cc=1;cc<=nmodes;cc++) {
    norm_tmp = (sqrt(sum(def(,,cc)^2)/norm))
    def(,,cc) /= norm_tmp; 
    gradG(,,cc,) /= norm_tmp;
  }
  
  return def;
}

func mygauss2(size,xc,yc,fwhmx,fwhmy,a,angle,fond,&grad,deriv=)
/* DOCUMENT mygauss2(size,xc,yc,fwhmx,fwhmy,a,angle,fond,&grad,deriv=)
* SEE ALSO:
*/
{
  x=(indgen(size))(,-:1:size);
  y=transpose(x);
  x -= xc;
  y -= yc;
  
  s = sin(angle);
  c = cos(angle);
  u = x*c-y*s;
  v = x*s+y*c;

  g1 = exp(-(u/(fwhmx/1.66))^2-(v/(fwhmy/1.66))^2);
  g = a*g1;
  
  if (deriv) {
    grad = array(0.,[3,size,size,7]);
    grad(,,1) = 2*(c*u/((fwhmx/1.66)^2)+s*v/((fwhmy/1.66)^2))*g; // gradxc
    grad(,,2) = 2*(c*v/((fwhmy/1.66)^2)-s*u/((fwhmx/1.66)^2))*g; // gradyc
    grad(,,3) = 5.5112*u^2*g/fwhmx^3; // grad fwhmx
    grad(,,4) = 5.5112*v^2*g/fwhmy^3; // grad fwhmy
    grad(,,5) = g1; // grad a
    grad(,,6) = 5.5112*g*u*v*(1./fwhmx^2-1./fwhmy^2); // grad angle
    grad(,,7) = 1; // grad fond
  }
  return g+fond;
}

func prepadonis(size, patchDiam, nactu, couplingFact, defpup, deftiptilt)
{
  pup = defpup;
  x=span(-size,size,size)(,-:1:size) / (patchDiam+1);
  y=transpose(x);
  x0 = span(-1,1,nactu)(,-:1:nactu);
  y0 = transpose(x0);
  k=0;
  espaceInterAct = 2./(nactu-1);
  def = array(0.0, size, size, nactu*nactu);
  A = -log(couplingFact)/(espaceInterAct^2);
  for(j=1;j<=nactu;j++) {
    for(i=1; i<=nactu;i++) {
      r = sqrt( x0(i,j)^2 + y0(i,j)^2 );
      if( r<1.1 ) {
        k++;
        def(,,k) = exp(-A*((x-1.*x0(nactu-i+1,j))^2+(y-1.*y0(i,j))^2));
        //def(,,k) = roll(def(,,k),[-1,-3]);
      }
    }
  }
  def = pup * def(,,:k);
  
  nmodes = dimsof(def)(4);
  norm = sum(deftiptilt(,,2)^2);
  for (cc=1;cc<=nmodes;cc++)
    def(,,cc) /= (sqrt(sum(def(,,cc)^2)/norm)); 
 
  return def;
}

func phase2psf(phase,mask,size,&ampliPup,&ampliFoc,fftws=)
/* DOCUMENT phase2psf
 * psf=phase2psf(phase,mask,size,&ampliPup,&ampliFoc,fftws=)
 * This routine returns the PSF estimated from a phase map
 *
 * KEYWORDS :
 * phase      (input)    : Phase map 
 * imwidth    (input)    : Number of pixels of the image
 * mask       (input)    : The pupil
 * ampliPup   (output)   : Complex amplitude in the pupil plane (size x size)
 *
 * SEE ALSO phase2psfMode
 */ 
{
  pupd = dimsof(mask)(2);
  sphase = dimsof(phase);

  if (sphase(2) != pupd) {
    write,"Size of Phase and Mask incompatible";
    error;
  }

  // Complex amplitude in the pupil plane
  ampliPup = array(complex,[2,size,size]);
  (ampliPup.re)(1:pupd,1:pupd) = cos(phase)*mask;
  (ampliPup.im)(1:pupd,1:pupd) = sin(phase)*mask;

  // Complex amplitude in the focal plane
  ampliFoc = fft(ampliPup,1,setup=fftws);
  psf      = abs(ampliFoc)^2;

  // normalization
  psf      = psf/sum(psf);

  return psf;
}


func pray_c_data(&grad_psf,&grad_obj,ft_object=,image=,ft_psf=,variance=,type=,nstar=,fftws=,grad_o=)
/* DOCUMENT pray_c_data
 * criterion=pray_c_data(&grad_psf,ft_object=,image=,ft_psf=,variance=,type=,fftws=)
 *
 * This routine returns the following criterion :
 * C(i) = sum (1/var x |(h * o) - i|^2)
 *
 * KEYWORDS :
 * grad_psf : The gradient of the criterion with respect to the psf
 * ft_object: The object FT (centered)
 * image    : The image 
 * ft_psf   : The psf FT (centered) 
 * variance : The noise variance. Can be a scalar or a 2D array
 *
 * SEE ALSO: 
 */ 
{
  extern pray_buffer,pray_ndefoc;
  extern pray_currerr;
  
  if (type == []) type = 0;
  if (grad_o == []) grad_o = 0;
  
  dim2 = numberof(image);
  np = long(sqrt(dim2));         //image must be square
  
  // Test on the variance
  if ((numberof(variance) != dim2) & (numberof(variance) != 1)) {
    error,"variance should be scalar or of same size than image.";
  }
  
  // estimation of (h*o-i)
  if (numberof(ft_object) > 1)
    HO = 1./dim2*(double(fft((ft_psf*ft_object),-1,setup=fftws)));
  else {
    HO = 1./dim2*roll(double(fft(ft_psf,-1,setup=fftws)));
    //renormalisation ... useful for point-like images
    norm_HO = sum(HO);
    HO /= norm_HO;
    HO *= sum(image);
  }
  
  HOminusI = HO - image;

  if (numberof(ft_object) > 1) 
    grad_psf = 1./dim2*double(fft(conj(ft_object)*fft(HOminusI/variance,1,setup=fftws),\
                                -1,setup=fftws));
  else grad_psf = 1./dim2*roll(double(fft(fft(HOminusI/variance,1,setup=fftws), \
                                          -1,setup=fftws)));
  
  if (grad_o)
    grad_obj = 1./dim2*double(fft(conj(ft_psf)*fft(HOminusI/variance,1,setup=fftws), \
                                       -1,setup=fftws));
  else grad_obj = 0.;
  
  if (type != []) {
    if (nstar != []) {
      pray_buffer(,,nstar,type+pray_ndefoc+1) = HO;
      pray_currerr(,,nstar,type) = (HOminusI)^2;
    } else {
      pray_buffer(,,type+pray_ndefoc+1) = HO;
      pray_currerr(,,type) = (HOminusI)^2;
    }
  }

  return .5 * sum((HOminusI)^2/variance);
}

func pray_gradpsf2param(gradPsf,&gradPhase,modesArray,mask,ampliPup,ampliFoc,pupd,fftws=,deltaFoc=,gradG=,coeffs=)
/* DOCUMENT 
 *  
 * grad=pray_gradpsf2param(gradPsf,&gradPhase,modesArray,mask,ampliPup,ampliFoc,pupd,fftws=,deltaFoc=)
 *
 * This routine returns the gradient with respect to the parameters (mode coeffs)
 * from the gradient with respect to the PSF
 *
 * KEYWORDS :
 * gradPsf   : The gradient with respect to the PSF
 * gradPhase : The gradient with respect to the Phase (pupSize x pupSize)
 * mask      : The pupil
 * modesArray: The modes (nbModes x pupsize^2)
 * ampliPup  : The complex amplitude in the pupil plane (pupSize x pupSize)
 *
 * SEE ALSO: 
 */ 
{
  // modesArray contains the Zi in column ...
  // in th tomographic case, modes Array contains the modes, as viewed in the pupil for
  // the specific direction ... need to use getModesInPup before
  //pupd2 = pud^2;
  size2 = numberof(mask);
  size  = sqrt(size2);
  
  dummy = where(mask != 0.);
  count = numberof(dummy);

  // here we compute the gradient of the psf with respect to the phase
  gradPhase = (-2./size2/count)*(ampliPup*fft(gradPsf*conj(ampliFoc),1,setup=fftws)).im;
  // gradphase is nil outside (3:pupd+2,3:pupd+2) so we need to recenter it for the
  // product with the modes
  gradPhase = roll(gradPhase,[size/2-pupd/2-2,size/2-pupd/2-2]);
  
  // here we compute the gradient of the psf with respect to the mode coeffs
  gradCoeff = modesArray(*,)(where(mask),)(+,)*gradPhase(*)(where(mask))(+);

  // here we compute the gradient of the psf with respect to the focus scale factor
  if (deltaFoc != []) gradScale = deltaFoc * ((modesArray(*,3)(where(mask)))(+)*\
                                              (gradPhase(*)(where(mask)))(+));

  // here we compute the gradient of the psf with respect to the gird parameters
  if (gradG != []) {
    gGrid = array(0.0,4);
    for (cc=1;cc<=4;cc++)     
      gGrid(cc) = (gradG(,,+,cc)*coeffs(+))(*)(where(mask))(+)* \
        gradPhase(*)(where(mask))(+);
  }
  
  return _(gGrid,gradScale,gradCoeff);
}

func pray_coeff2psfs(coeff,&ampliPup,&ampliFoc,poffset=)
/* DOCUMENT pray_coeff2psfs
 *  
 * psfs = pray_coeff2psfs(coeff,&ampliPup,&ampliFoc,scale=)
 *
 * This routine returns a PSF array from an array of coefficients
 *
 * KEYWORDS :
 * 
 * coeff   : The coefficients
 * ampliPup: The complex amplitude in the pupil plane 
 * ampliFoc: The complex amplitude in the focal plane 
 * scale   : an optional scale on the first coefficient
 * poffset : a phase map or a scalar used as an offset to the phase
 *
 * SEE ALSO: 
 */ 
{
  extern pray_data,pray_mircube;

  mircube     = *(pray_data.mircube);
  def         = *(pray_data.def);
  ipupil      = *(pray_data.ipupil);  
  xshifts     = *(pray_data.xshifts);
  yshifts     = *(pray_data.yshifts);
  ishifts     = *(pray_data.ishifts);
  jshifts     = *(pray_data.jshifts);
  _n          = pray_data._n;
  _n1         = pray_data._n1;
  _n2         = pray_data._n2;
  nzer        = *(pray_data.nzer);
  alt         = *(pray_data.alt);
  size        = pray_data.size;
  fftws       = *(pray_data.fft_ws);
  fit_starpos = pray_data.fit_starpos;
  
  nalt     = numberof(alt);
  cent     = size/2+0.5;
  ntarget  = dimsof(xshifts)(4);
  
  mircube *= 0.0f;
  istart = 0;

  if (fit_starpos) {
    coeff_starpos = coeff(1:2*ntarget); // organised as xy xy xy ....
    coeff = coeff(2*ntarget+1:);
  }
  
  for (i=1;i<=nalt;i++) {
    if ((fit_starpos) && (i>1)) tmp = coeff(istart-1:nzer(i)+istart-2);
    else  tmp = coeff(istart+1:nzer(i)+istart);
    if (i==1) {
      if (fit_starpos) {
        tmp2 = def(,,istart+1:nzer(i)+istart)(,,3:);
        tmp = coeff(istart+1:nzer(i)-2+istart);
      } else tmp2 = def(,,istart+1:nzer(i)+istart);
    } else tmp2 = def(,,istart+1:nzer(i)+istart);
    mircube(,,i) = float(tmp(+)*tmp2(,,+));
    istart += nzer(i);
  }
  
  pray_data.mircube = &float(mircube);
  
  bphase = array(float,[3,size,size,ntarget]);
  psnx = dimsof(mircube)(2);
  psny = dimsof(mircube)(3);
  
  for (k=1;k<=ntarget;k++) {
    sphase = array(float,_n,_n);  
    if (fit_starpos) {
      mircube2 = mircube;
      scoeffs = coeff_starpos(2*k-1:2*k);
      mircube2(,,1) += scoeffs(+)*def(,,1:2)(,,+)
        err = _get2dPhase(&float(mircube2),psnx,psny,nalt,&sphase,_n,_n,&int(ishifts(,,k)), \
                  &float(xshifts(,,k)),&int(jshifts(,,k)),&float(yshifts(,,k)));
     } else err = _get2dPhase(&mircube,psnx,psny,nalt,&sphase,_n,_n,&int(ishifts(,,k)), \
                      &float(xshifts(,,k)),&int(jshifts(,,k)),&float(yshifts(,,k)));
    if (err != 0) {error,"Error in getPhase2dFromDms";}
    bphase(_n1:_n2,_n1:_n2,k) = float(sphase);
  }
  
  psfs = array(float,[3,size,size,ntarget]);
  ampliPup = ampliFoc = array(complex,[3,size,size,ntarget]);

  for (i=1;i<=ntarget;i++) {
    tmp = bphase(,,i);
    if (poffset != []) tmp -= poffset;
    psfs(,,i) = phase2psf(tmp(_n1:_n2,_n1:_n2),ipupil(_n1:_n2,_n1:_n2),\
                          size,amp1,amp2,fftws=fftws);
    ampliPup(,,i) = amp1;
    ampliFoc(,,i) = amp2;
  }
  return psfs;
}

func pray_init(xpos,ypos,tmodes=,tiptilt=,poffset=)
/* DOCUMENT pray_init
 *  
 * pray_init,xpos,ypos,tmodes=,tiptilt=
 *
 * In this routine, the overall geometry is initialized
 *
 * KEYWORDS :
 * 
 * xpos   : x positions of the various targets
 * ypos   : y positions of the various targets
 * tmodes : the kind of modes (zernike, kl, mirror)
 * tiptilt: a flag to include or not tip-tilt
 *
 * SEE ALSO: 
 */ 
{
  extern pray_data,pray_mircube,def;

  teldiam = pray_data.teldiam;
  cobs    = pray_data.cobs;
  pupd    = pray_data.pupd;
  size    = pray_data.size;
  alt     = *(pray_data.alt);
  nzer    = *(pray_data.nzer);
  
  cent    = size/2+0.5;
  nalt    = numberof(alt);
  ntarget = numberof(xpos);

  // Create mircube
  mircube = array(float,[3,size,size,nalt]);
  pray_data.mircube = &mircube;

  _p  = pupd;
  _p1 = long(ceil(cent-_p/2.));
  _p2 = long(floor(cent+_p/2.));
  _p  = _p2-_p1+1;
    
  _n  = _p+4;   pray_data._n = _n;
  _n1 = _p1-2;  pray_data._n1 = _n1;
  _n2 = _p2+2;  pray_data._n2 = _n2;

  // Init ipupil
  ipupil = float(make_pupil(size,pupd,xc=cent,yc=cent,cobs=cobs));
  pray_data.ipupil = &ipupil;

  // Init def
  def = array(float,[3,size,size,nzer(sum)]);
  
  psize = teldiam/pupd;
  cpt = 0;

  nact_gems = [17,20,12];
  
  if (tmodes == 5) {
    nxact_dm = [19,22,16];
    pitch_dm = [float(pupd)/nxact_dm(1),float(pupd)/nxact_dm(1),2.*pupd/nxact_dm(1)];
    patchDiam_dm = pitch_dm * nxact_dm;
    coupl = [0.33,0.33,0.2];
    pitchMarg = [0.5,1.2,0.5];
    unitV = [0.6*1.308,0.6*1.224,0.6*1.123];
  }
    
  // boucle sur les couches
  for (k=1;k<=nalt;k++) {
    patchDiam = long(pupd+2.*max(abs(xpos,ypos))*4.848e-6*(alt(k))/psize);
    if ((patchDiam % 2) != 0) patchDiam += 1;
    // taille de la pupille sur le support
    prepzernike,size,patchDiam,cent,cent; // size = taill tot support, cent=size/2+0.5
    if( tmodes == 0 ) {   
      if (k == 1) {
        //selz = _(4,2,3,4+indgen(nzer(k)));
        for (i=1;i<=nzer(k);i++) {
          cpt ++;
          def(,,cpt) = zernike(i+1);
          //def(,,cpt) = zernike(selz(i));
        }
      } else  {
        for (i=7;i<=nzer(k)+6;i++) {
          cpt ++;
          def(,,cpt) = zernike(i);
        }
      }
    }

    if(tmodes == 1) {
      pup1 = [];
      if (k==1) {
        foc = zernike(4);
        kl = make_kl(nzer(k),patchDiam,v,obas,pup1,oc=0.0,nr=128);
        //    kl = order_kls(kl,patchDiam,upto=20);
        def(size/2-patchDiam/2+1:size/2+patchDiam/2,            \
            size/2-patchDiam/2+1:size/2+patchDiam/2,2:nzer(k)) = kl(,,:-1);
        
        //def = def(,,_(1,indgen(nzer(k)-1)));
        def(,,1:nzer(k)) *= ipupil;
        //     m = def(*,)(+,) * foc(*)(+);   
        //     def -= m(+)*foc(,,-:1:nzer(k))(,,+);
        def(,,1) = foc;
        cpt += nzer(k);
      } else {
        if ((patchDiam % 2) != 0) patchDiam += 1;
        kl = make_kl(nzer(k)+5,patchDiam,v,obas,pup1,oc=0.0,nr=128);
        def(size/2-patchDiam/2+1:size/2+patchDiam/2,            \
            size/2-patchDiam/2+1:size/2+patchDiam/2,cpt+1:nzer(k)+cpt) = kl(,,6:);
        def(,,cpt+1:nzer(k)+cpt) *= ipupil;
        cpt += nzer(k);
      }
    } 
    
    if(tmodes == 2) {
      // rough mirror model
      deftt = [zernike(4),zernike(2),zernike(3)];
      tmp = prepadonis(size, pupd, nact_gems(k), 0.2, ipupil, deftt);
      if (k==1) {
        def(,,1:3) = [zernike(2),zernike(3),zernike(4)];
        def(,,4:nzer(k)) = tmp;
        istart = nzer(k)
      } else {
        def(,,istart+1:nzer(k)+istart) = tmp;
        istart += nzer(k);
      }
      // size = taill tot support, cent=size/2+0.5
    }
    
    if(tmodes == 3) {
      nactu = ceil(sqrt(nzer(1)));
      espaceInterAct = 2./(nactu-1);
      dmpitch = espaceInterAct*pupd/2.;
      x0 = (size-pupd+1)/2.-dmpitch; // this seems like a normal initial value
      y0 = (size-pupd+1)/2.-dmpitch; 
      // grid model
      def(,,1) = zernike(4);
      def(,,2) = zernike(2);
      def(,,3) = zernike(3);
      def(,,4:) = prepmirror(size,pupd,nactu,0.09,dmpitch,dmpitch,x0,y0,ipupil,def(,,1:3),&gradG,orthonorm=1);
    }

    // 1 piston / 2 tt / 3 quadratic / defoc = 5
    if( tmodes == 4 ) {
      if (k == 1) { // include tt
        istart = 0;
        tmp_def = make_diskharmonic(size,patchDiam,nzer(k)+1,xc=cent,yc=cent)(,,2:);
        tmp = tmp_def(,,3); // putting back focus where it belongs
        tmp_def(,,3) = tmp_def(,,4);
        tmp_def(,,4) = tmp;
        tmp = [];
        def(,,istart+1:nzer(k)+istart) = tmp_def;
      } else  {
       tmp_def = make_diskharmonic(size,patchDiam,nzer(k)+6,xc=cent,yc=cent)(,,7:);
       def(,,istart+1:nzer(k)+istart) = tmp_def;
      }
      istart += nzer(k);
    }
    if( tmodes == 5 ) { // using GeMS DMs
      //make_pzt_dm(size,dim,nxact,cobs,pitch,coupling,xflip=,yflip=,pitchMargin=,unitpervolt=)
      tmp = make_pzt_dm(size,patchDiam_dm(k),nxact_dm(k),pray_data.cobs,pitch_dm(k),
                        coupl(k),pitchMargin=pitchMarg(k),unitpervolt=unitV(k));
      
      if (k==1) {
        prepzernike,size,patchDiam_dm(k),cent,cent; 
        def(,,1:3) = [zernike(2),zernike(3),zernike(4)];
        tmp2 = filter_quad(size,patchDiam,cent,cent,tmp,foconly=1);
        def(,,4:nzer(k)) = tmp2;
        istart = nzer(k);
      } else {
        tmp2 = filter_quad(size,patchDiam,cent,cent,tmp);
        def(,,istart+1:nzer(k)+istart) = tmp;
        istart += nzer(k);
      }
    }
  }

  pray_data.def = &def;
  pray_mircube = mircube*0.0f;

  // Init dmgsXYposcub : useful for get2dPhase
  xref = indgen(_n)-(_n+1)/2.;
  yref = indgen(_n)-(_n+1)/2.;
  
  dmgsxposcub = dmgsyposcub = xshifts = yshifts = ishifts = jshifts = \
    array(float,[3,_n,nalt,ntarget]);
  
  for (n=1;n<=ntarget;n++) {
    // loop on pseudo-DMs
    for (ns=1;ns<=nalt;ns++) {
      // offsets of the center of beam on DM NS
      xoff = xpos(n)*4.848e-6*alt(ns)/psize;
      yoff = ypos(n)*4.848e-6*alt(ns)/psize; 
      dmgsxposcub(,ns,n) = xref + xoff;
      dmgsyposcub(,ns,n) = yref + yoff;
    }   
  }

  for (k=1;k<=ntarget;k++) {
    xshifts(,,k) = dmgsxposcub(,,k)+(cent-1)(-,);
    yshifts(,,k) = dmgsyposcub(,,k)+(cent-1)(-,);  
    ishifts(,,k) = int(xshifts(,,k)); xshifts(,,k) = xshifts(,,k) - ishifts(,,k);
    jshifts(,,k) = int(yshifts(,,k)); yshifts(,,k) = yshifts(,,k) - jshifts(,,k);
  }
  
  pray_data.ishifts = &ishifts;
  pray_data.jshifts = &jshifts;
  pray_data.xshifts = &xshifts;
  pray_data.yshifts = &yshifts;

  szdef = dimsof(def);
  defPup = array(float,[4,szdef(2),szdef(3),szdef(4),ntarget]);

  for (i=1;i<=ntarget;i++) defPup(,,,i) = getDefInPupilFromDir(i);

  pray_data.defPup = &defPup; 
}


func getDefInPupilFromDir(ndir)
/*
  return the array of modes as viewed from the pupil in a specified direction
 */
{
  extern pray_data;

  mircube     = *(pray_data.mircube);
  alt         = *(pray_data.alt);
  nzer        = *(pray_data.nzer);
  def         = *(pray_data.def);
  xshifts      = *(pray_data.xshifts);
  yshifts      = *(pray_data.yshifts);
  ishifts      = *(pray_data.ishifts);
  jshifts      = *(pray_data.jshifts);
  size        = pray_data.size;
  _n          = pray_data._n;
  _n1         = pray_data._n1;
  _n2         = pray_data._n2;

  cent        = size/2+0.5;
  nalt        = numberof(alt);

  psnx        = dimsof(mircube)(2);
  psny        = dimsof(mircube)(3);

  // geometry init : get the proper points coordinates : dmgsxposcub
  mydef = def*0.;

  mircube *= 0.0f;
  cpt      = 0; 
  for (i=1;i<=nalt;i++) {
    for (j=1;j<=nzer(i);j++) {
      mircube *= 0.0f;
      bphase = array(float,[2,size,size]);
      sphase = array(float,_n,_n);
      cpt++;
      mircube(,,i) = float(def(,,cpt));
      err =                                                             \
        _get2dPhase(&mircube,psnx,psny,nalt,&sphase,_n,_n,&int(ishifts(,,ndir)),&float(xshifts(,,ndir)),&int(jshifts(,,ndir)),&float(yshifts(,,ndir)));
      if (err != 0) {error,"Error in getPhase2dFromDms";}
      bphase(_n1:_n2,_n1:_n2) = float(sphase);
      mydef(,,cpt) = bphase;
    }
  }
  return mydef;
}

func pray_error(param,&gradient,extra)
/* DOCUMENT pray_error
 *  
 * criterion=pray_error(param,gradient,extra)
 *
 * This routine returns the error to be minimized in pray
 *
 * KEYWORDS :
 * param   : The parameters to be mimimized
 * gradient: The gradient of this error with repect to the parameters
 * extra   : A structure containing all necessary parameters
 *           (see the content of pday_struct() for more details).
 *
 * SEE ALSO: 
 */ 
{ 
  extern pray_mircube;
  local deltaFoc;
  
  images      = *extra.images;
  variance    = *extra.variance;
  ftobject    = *extra.ftobject;
  norm        = *extra.norm;
  deltaFoc    = *extra.deltaFoc;
  nzer        = *extra.nzer;
  alt         = *extra.alt;
  def         = *extra.def;
  defPup      = *extra.defPup;
  ipupil      = *extra.ipupil;
  pupd        = extra.pupd;
  size        = extra.size;
  fftws       = *extra.fft_ws;
  scale       = extra.scale;
  dmgrid      = extra.dmgrid;
  poffset     = *extra.poffset;
  diff_tt     = extra.diff_tt;
  fit_starpos = extra.fit_starpos
  fit_object  = extra.fit_object

  if (dmgrid) {
    // here we need to recompute the modeArray + gradients
    def      = *(pray_data.def);
    x0       = param(1);
    y0       = param(2);
    dmpitchX = param(3);
    dmpitchY = param(4);
    gradG = 0.0;
    def(,,4:) = prepmirror(size,pupd,ceil(sqrt(nzer(1))),0.09,dmpitchX,dmpitchY,x0,y0, \
                           ipupil,def(,,1:3),gradG,orthonorm=1) ;
    param         = param(5:);
    pray_data.def = &float(def);
  } else gradG = [];

  if (fit_object) {
    params_object = param(1:3);
    param = param(4:);
  }
  // we deal with scale & diff_tt here
  // we deal with starpos in coeff2psfs (faster)
  if (scale) {
    scale = param(1);
    coeffs   = param(2:);
  } else {
    scale = [];
    coeffs   = param;
  }

  if (diff_tt) {
    coeff_diff = coeffs(1:2*numberof(deltaFoc)); // organized as xy xy xy ...
    coeffs = coeffs(2*numberof(deltaFoc)+1:);
  }
  
  ntarget = dimsof(images)(4);

  nModes = nzer(*)(sum);

  // int crit_array
  crit_array = array(0.0f,(numberof(deltaFoc)+1)*ntarget);
  
  // scale gradient has to be summed on all defoc images and all stars (same coeff for everyone)
  // diff_tt gradient has to be summed on all stars for a specific defoc pos
  // starpos gradient has to be summed on all defoc pos of each star
  
  if (fit_starpos) gradientModes = array(0.0,nModes-2); // removing global tt
  else gradientModes = array(0.0,nModes);
  if (scale != []) gradientScale = 0.0f;
  if (dmgrid) gradientGrid = array(0.0,4,(numberof(deltaFoc)+1)*ntarget);
  if (fit_starpos) gradientStars = array(0.0,2,ntarget);
  if (diff_tt) gradientDiff = array(0.0,2,numberof(deltaFoc));
  
  if (fit_object) {
    gradientObj = array(0.0,3);
    new_obj = mygauss2(size,size/2+1,size/2+1,params_object(1),params_object(2),1.,params_object(3),0.,grad_obj,deriv=1);
    norm_obj = sum(new_obj);
    new_obj /= norm_obj;
    new_obj *= images(,,avg,1)(*)(sum);
    grad_obj /= norm_obj;
    grad_obj *= (images(,,avg,1)(*)(sum));
    ftobject = fft(new_obj,1);
  }
  
  //-----------------------------------------------------------------
  // First dealing with in-focus images
  tmp = coeffs;
  psfs = pray_coeff2psfs(tmp,ampliPup,ampliFoc,poffset=poffset);
  // here the psf is centered at (0,0)
  // tmp includes tt coeffs for starpos estimation

  for (i=1;i<=ntarget;i++) {
    ftPsf = fft(psfs(,,i),1,setup=fftws);
    // Estimation of the criterion associated with image #i
    if (numberof(variance) == 1) var = variance;
    else {
      if (numberof(variance) == numberof(deltaFoc)+1) var = variance(i);
      else var = variance(,,i,1);
    }

    nstar = (ntarget > 1 ? i : []);
    crit_array(i) =
      pray_c_data(gradientPsf,gradObj,ft_object=ftobject,image=images(,,i,1),   \
                  ft_psf=ftPsf,variance=var,type=1,nstar=nstar,fftws=fftws,grad_o=fit_object);
    
    delta = []; // here there is no scale !
    
    tmp = pray_gradpsf2param(gradientPsf,gradPhaseOut,defPup(,,,i),ipupil,\
                             ampliPup(,,i),ampliFoc(,,i),pupd,fftws=fftws,\
                             deltaFoc=delta,gradG=gradG,coeffs=coeffs(4:));
    if (dmgrid) {
      gradientGrid(,i) = tmp(1:4);
      tmp = tmp(5:);
    }
    
    if (fit_starpos) gradientModes += tmp(3:);
    else gradientModes += tmp;
    
    if (fit_starpos) gradientStars(,i) = tmp(1:2);

    if (fit_object) {
      gradientObj(1:2) += grad_obj(*,3:4)(+,)*gradObj(*)(+);
      gradientObj(3) += grad_obj(*,6)(+)*gradObj(*)(+);
    }
  }
  
  pray_mircube = (*pray_data.mircube)*ipupil;

  //----------------------------------------------------------------
  //Extra-Focal images (we can introduce as many images as we want)
  for(n=1;n<=numberof(deltaFoc);n++) {
    tmp = coeffs;
    if (scale != []) {
      if (fit_starpos) tmp(2*ntarget+1) += deltaFoc(n)*scale;
      else tmp(3) += deltaFoc(n)*scale;
    } else {
      if (fit_starpos) tmp(2*ntarget+1) += deltaFoc(n);
      else tmp(3) += deltaFoc(n);
    }

    
    if (diff_tt) {
      tmp(1) += coeff_diff(2*n-1);
      tmp(2) += coeff_diff(2*n);
    }
       
    psfs = pray_coeff2psfs(tmp,ampliPup,ampliFoc,poffset=poffset);
 
    for (i=1;i<=ntarget;i++) {
      ftPsf = fft(psfs(,,i),1,setup=fftws);
      // Estimation of the criterion associated with image #i
      if (numberof(variance) == 1) var = variance;
      else {
        if (numberof(variance) == numberof(deltaFoc)+1) var = variance(i);
        else var = variance(,,i,n+1);
      }
      nstar = (ntarget > 1 ? i : []);
      crit_array(n*ntarget+i) = pray_c_data(gradientPsf,gradObj,ft_object=ftobject, \
                                            image=images(,,i,n+1),ft_psf=ftPsf, \
                                            variance=var,type=1+n,nstar=nstar,fftws=fftws,grad_o=fit_object);
      if (scale != []) delta = deltaFoc(n);
      else delta = [];
      
      tmp = pray_gradpsf2param(gradientPsf,gradPhaseOut,defPup(,,,i),ipupil,\
                               ampliPup(,,i),ampliFoc(,,i),pupd,fftws=fftws,\
                               deltaFoc=delta,gradG=gradG,coeffs=coeffs(4:));
      if (dmgrid) {
        gradientGrid(,i+ntarget*n) = tmp(1:4);
        tmp = tmp(5:);
      }
      if (scale != []) {
        gradientScale += tmp(1);
        tmp = tmp(2:);
      }
      
      if (fit_starpos) gradientModes += tmp(3:);
      else gradientModes += tmp;
    
      if (fit_starpos) gradientStars(,i) += tmp(1:2);
      if (diff_tt) gradientDiff(,n) += tmp(1:2);
      
      if (fit_object) {
        gradientObj(1:2) += grad_obj(*,3:4)(+,)*gradObj(*)(+);
        gradientObj(3) += grad_obj(*,6)(+)*gradObj(*)(+);
      }
    }
  }
  //----------------------------------------------------------------

  gradientModes = gradientModes(*);
  if (fit_starpos) gradientStars = gradientStars(*);
  if (diff_tt) gradientDiff = gradientDiff(*);
  if (dmgrid) gradientGrid = gradientGrid(,sum);
  gradient = double(_(gradientGrid,gradientObj,gradientScale,gradientDiff,gradientStars,gradientModes));
  //gradient*=norm;
  
  return double(crit_array(sum));
}
 

func pray(images,xpos,ypos,deltaFoc,sigma,object,lambda,nzer=,alt=,teldiam=,cobs=,osampl=,disp=,verbose=,threshold=,nbiter=,tmodes=,tiptilt=,variance=,guess=,scale=,diff_tt=,fit_starpos=,fit_object=)
/* DOCUMENT pray
 *  
 *  pray(images,xpos,ypos,deltaFoc,sigma,object,lambda,nzer=,alt=,teldiam=,cobs=,osampl=,disp=,verbose=,threshold=,nbiter=,tmodes=,tiptilt=,variance=,guess=,scale=)
 *
 * KEYWORDS :
 * images  : The images (a 2D arrays) 
 * guess   : A guess
 * cobs    : The size of central obscuration on the pupil
 * sigma   : The noise standard deviation in the image
 * nbiter  : Number max of iterations for the criterion minimization
 * disp    : A flag to display (or not) the step-by-step results
 * verbose : A flag to print (or not) the details 
 * SEE ALSO:
 */ 
{
  extern pray_data,pray_selected_display;
  extern pray_iter,pray_param,dispok;
  extern stop_pray_loop;
  extern elapsed,cpu_start,pray_task,pray_eval,pray_step,pray_ws;
  
  local variance;
  
  // Init verbose and disp flags
  if (!is_set(verbose)) verbose=0;
  if (!is_set(disp)) disp=0;                           \
                               
  // sizes init and check
  dims = dimsof(images);
  if (dims(1) != 4) {
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","images format must follow : [size,size,ntarget,ndefoc]");
    return 1;
  }
  size = dims(2);
  ntarget = dims(4);
  ndefoc = dims(5)-1;
  if ((numberof(xpos) != ntarget) || (numberof(ypos) != ntarget)) {
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","incompatible images and xpos size");
    return 1;
  }

  // geometry init
  if (teldiam == []) teldiam = 7.9;
  if (cobs == [])    cobs    = 0.12;
  //osampl = (2.16e-6/8./4.848e-6)/pixSize;
  if (osampl == []) osampl = 1.;
  pupd=long(floor(size/2/osampl));
  if (pupd%2 != 0) pupd += 1;

  pray_data         = pray_struct();
  pray_data.teldiam = teldiam;
  pray_data.cobs    = cobs;
  pray_data.pupd    = pupd;
  pray_data.size    = size;

  if (poffset != []) {
    if ((numberof(poffset)!=size*size) || (numberof(poffset)!=1))
      error,"wrong size for poffset";
    else pray_data.poffset = &poffset;
  }

  // init nzer and alt if needed
  if (nzer == []) nzer = [10,10,10];
  if (alt == [])  alt  = [0.,4500.,9000.];

  if (tmodes == []) tmodes = 0;
  
  if (tmodes==3) dmgrid = 1;
  else dmgrid = 0;

  pray_data.nzer    = &nzer;
  pray_data.alt     = &alt;
  pray_data.dmgrid  = dmgrid;

  // pray init
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",\
             "Initializing the pray workspace ...");
  
  pray_init,xpos,ypos,tmodes=tmodes;
  
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Done ...");

  // ... testing the size and content of variance
  if (variance==[]) {
    variance = sigma^2;
    sz_variance = 1;
  } else {
    sz_variance = numberof(variance);
    if ((sz_variance != 1) & (sz_variance != numberof(deltaFoc)) & (sz_variance != numberof(images))) {
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",\
                 "variance must be a scalar or of same size than image.");
      return 1;
    } else {
      if (min(variance) <= 0.) {
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",\
                   "variance must be strictly positive");
        return 1;
      }
    }
  }
  
  // ... testing guess and create one if nill
  norm = 0.;//guess*0.;
  if (guess != []) {
    if (numberof(guess) != nzer(*)(sum)) {
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",\
                 "Guess does not have the proper size. Initializing to 0s.");
      pray_param = array(0.,nzer(*)(sum)) ;
    } else pray_param = guess;
  } else pray_param = array(0.,nzer(*)(sum)) ;
  //pray_param = (gaussdev(nzer(sum))*1./sqrt(indgen(nzer(sum))));

  if (!is_set(scale)) scale=0;
  if (!is_set(fit_starpos)) fit_starpos =0;
  if (!is_set(diff_tt)) diff_tt=0;
  if (!is_set(fit_object)) fit_object=0;

  if (fit_starpos) pray_param = _(array(0.0f,2*ntarget-2),pray_param);
  //-2 because we don't take into account the global tt
  if (diff_tt) pray_param = _(array(0.0f,2*ndefoc),pray_param);
  if (scale) pray_param = _(1.,pray_param);
  if (fit_object) pray_param = _([2.,2.,0.],pray_param);

  
  // adding parameters for the dm grid model
  if (dmgrid) {
    nactu = ceil(sqrt(nzer(1)));
    espaceInterAct = 2./(nactu-1);
    dmpitch = espaceInterAct*pupd/2;
    x0 = (size-pupd+1)/2.-dmpitch; // this seems like a normal initial value
    y0 = (size-pupd+1)/2.-dmpitch;
    pray_param = _(x0,y0,dmpitch,dmpitch,pray_param);
  }

  if( numberof(object) == 1 ) ftobject = 1.0f;
  else {
    if (fit_object) object = mygauss2(size,size/2+1,size/2+1,1.0,1.0,1.,0,0.);
    ftobject=fft(object,1);
  }

  if (numberof(threshold) == 0) {
    if (typeof(images) == "float") {
      mythresh = 1e-7;
    } else mythresh = 1e-16; 
  } else mythresh = threshold;


  pray_data.ftobject    = &ftobject;
  pray_data.variance    = &variance;
  pray_data.images      = &images;
  pray_data.norm        = &norm;
  pray_data.deltaFoc    = &deltaFoc;
  //pray_data.fft_ws     = &(fft_setup(dimsof([size,size])));
  pray_data.threshold   = mythresh;
  pray_data.lambda      = lambda;
  pray_data.nbiter      = nbiter;
  pray_data.neval       = 10*nbiter;
  pray_data.scale       = scale;
  pray_data.fit_starpos = fit_starpos;
  pray_data.diff_tt     = diff_tt;
  pray_data.fit_object  = fit_object;

  /*
  // Some initializations for the minimization process
  for (cc=numberof(deltaFoc)+1;cc<=numberof(deltaFoc)+dimsof(pray_mircube)(4);cc++) {
    pname = swrite(format="Phase on layer %d",cc-numberof(deltaFoc));
    pyk,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",\
               cc,pname);
  }
  pyk,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_active(%d)",0);
  pray_selected_error=1;
  */
  pray_update_zernike_table,pray_param,pray_data.lambda,nzer(1);
  pyk,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_active(%d)",0);
  pray_selected_display=1;

  pray_param = double(pray_param); //optimPack compat ... needs to be double
  nmin = numberof(pray_param);
  min_method = 0;
  
  if (min_method == 0) {
    /* Variable metric. */
    if (is_void(ndir)) ndir = nmin;
    method_name = swrite(format="Limited Memory BFGS (VMLM with NDIR=%d)",
                         ndir);
    pray_ws = op_vmlmb_setup(nmin, ndir, fmin=0.0);
  } else if (min_method < 0) {
    if (is_void(ndir)) ndir = nmin;
    method_name = swrite(format="Limited Memory BFGS (LBFGS with NDIR=%d)",
                         ndir);
    pray_ws = op_lbfgsb_setup(nmin, ndir);
  } else if (min_method >= 1 && min_method <= 15) {
    /* Conjugate gradient. */
    ndir = 2;
    method_name = swrite(format="Conjugate Gradient (%s)",
                         ["Fletcher-Reeves", "Polak-Ribiere",
                          "Polak-Ribiere with non-negative BETA"](min_method&3));
    pray_ws = optim_cgmn_setup(min_method);
  }
  
  stop_pray_loop = 0;
  dispok = disp;
  
  pray_step = 0.0;
  pray_task = 1;
  pray_eval = pray_iter = 0;
  elapsed = array(double, 3);
  timer, elapsed;
  cpu_start = elapsed(1);

  pray_loop;
}

func pray_loop(one)
{
  extern pray_data,pray_iter;
  extern pray_param,dispok;
  extern stop_pray_loop;
  extern elapsed,cpu_start,pray_task,pray_eval,pray_step,pray_wfs,pray_f,pray_g;

  if (pray_task == 1) {
    /* Evaluate function. */
    pray_f = pray_error(pray_param,pray_g,pray_data);
    ++pray_eval;
  }
  
  pray_iter+=1;
  /* Check for convergence. */
  if (pray_task != 1 || pray_eval == 1) {
    timer, elapsed;
    cpu = elapsed(1) - cpu_start;
    gnorm = sqrt(sum(pray_g*pray_g));
    
    pray_update_comments,pray_iter, pray_eval, cpu,pray_f,gnorm,pray_step;
    
    if (pray_iter == pray_data.nbiter) {
      pray_progressbar_text,"I reached the number max of iterations";
      pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                 "End of The Minimization process");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",          \
                 "---------------------------------------------------------------");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",  \
                 "I reached the number max of iterations");
      stop_pray_loop = 1;
    }

    if (pray_eval == pray_data.neval) {
      pray_progressbar_text,"I reached the number max of eavaluation";
      pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                 "End of The Minimization process");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",          \
                 "---------------------------------------------------------------");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",  \
                 "I reached the number max of evaluations");
      stop_pray_loop = 1;
    }
    if (dispok) {
      pray_disp_select;
      pray_disp_error;
    }

    pray_progressbar_frac,float(pray_iter)/pray_data.nbiter;
    pray_progressbar_text,swrite(format="iter %d of %d",pray_iter,pray_data.nbiter);
  }

  min_method = 0;
  
  if (min_method == 0) {
      pray_task = op_vmlmb_next(pray_param, pray_f, pray_g, pray_ws);
      pray_iter = (*pray_ws(2))(7);
      pray_step = (*pray_ws(3))(22);
  } else {
    if (min_method < 0) {
      pray_task = op_lbfgsb_next(double(pray_param), pray_f, pray_g, pray_ws);
      if (pray_task == 2 || pray_task == 3) ++pray_iter;
      pray_step = -1.0;
    } else {
      optim_cgmn, pray_param, pray_f, pray_g, pray_task, pray_ws;
      pray_iter = optim_cgmn_iter(pray_ws);
      pray_step = optim_cgmn_step(pray_ws);
    }
  }
  
  if (pray_task > 2) {
    pray_progressbar_frac,1.0;
    pray_progressbar_text,"I reached the convergence threshold";
    pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')","End of The Minimization process");
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","I reached the convergence threshold");
    stop_pray_loop = 1;
  }    
 
  pray_update_zernike_table,pray_param,pray_data.lambda,(*pray_data.nzer)(1);

  if (stop_pray_loop) {
    stop_pray_loop=0;
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
    newstring = swrite(format="Value of the criterion after %d iterations : %23.2e",pray_iter,pray_error(pray_param,grad_fin,pray_data));
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
    
    pray_disp_select;
    pray_disp_error;
    pyk,cmd_pray+"y_on_loop_finished()";
    return;
  }
  
  if (!one) {
    after,0.1,pray_loop;
  }
}

func pray_update_comments(iter, eval, cpu,f,gnorm,step)
{
  extern cmd_pray;
  
  pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')"," ITER    EVAL     CPU [s]          FUNC          GNORM   STEPLEN  ");
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","------    ------     ----------   -----------------------    ---------     ---------");
  newstring = swrite(format=" %5d %5d %10.3f     %+-10.6e  %-9.1e  %-9.1e",iter, eval, cpu,f , gnorm, step);
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
}

func pray_print_results(output, iter, eval, cpu, fx, gnorm, steplen, x, extra)
{
  extern cmd_pray;
  
  pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')","EVAL CPU[ms]               FUNC                     GNORM   STEPLEN");
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","--------------------------------------------------------------------------------");
  newstring = swrite(format="%4d %10.1f  %+-24.15e%-9.1e%-9.1e",eval, cpu, fx, gnorm, step);
  pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
}

func pray_update_zernike_table(coeffs,lambda,nzer)
{
  //box_content = swrite(format="%4.1f",deltaFoc/(2*pi/lambda));
  //pyk,swrite(format=cmd_pray+"glade.get_widget('deltafoc_entry').set_text('%s')",box_content);ic
  istart = 1;
  if (pray_data != []) {
    if (pray_data.dmgrid) coeffs = coeffs(5:);
    if (pray_data.scale) coeffs = coeffs(2:);
    if (pray_data.diff_tt) coeffs = coeffs(2*numberof(*pray_data.deltaFoc)+1:);
    if (pray_data.fit_starpos) {
      coeffs = coeffs(2*dimsof(*pray_data.images)(4)+1:);
      istart = 3;
    }
    if (max(coeffs) == 0.) istart = 1;
  }
  
  for (i=istart;i<=clip(nzer,,20);i++) {
    box_name = swrite(format="mod%d",i);
    box_content = swrite(format="%4.1f",coeffs(i+1-istart)/(2*pi/lambda));
    pyk,swrite(format=cmd_pray+"glade.get_widget('"+box_name+"').set_text('%s')",box_content);
  } 
}


 /*
   
   system oriented phase diversity

   first we determine a mirror model from in & out of focus images in close loop,
   with some waffle on the dm. waffle should help us to get a good accuracy to measure
   the 4 dm model parameters : x0, y0, dmpitchX and dmpitchY

   then we acquire one set of phase diversity images per actuator with each time
   one actuator pushed. we reconstruct the set of coefficients. this vector when
   clipped is one for the pushed actuator and 0 everywhere else. using all the
   phase diversity sets we build a matrix that maps the real dm actuators onto the
   simulated dm actuators.

   from that point, we can easily iterate with classical phase diversity using the
   above determined dm model

  */

func pray_script_mirror(targetx,targety,nactu,deltaFoc_nm,lambda_im,teldiam,cobs,pix_size,obj_type,obj_size,tiptilt,size,snr)
{
  
  // create psfs
  espaceInterAct = 2./(nactu-1);
  dmpitch = espaceInterAct*pupd/2;
  x0 = (size-pupd+1)/2.; // this seems like a normal initial value
  y0 = (size-pupd+1)/2.; 

  cent    = size/2+0.5;
  ipupil = float(make_pupil(size,pupd,xc=cent,yc=cent,cobs=cobs));
  def	= array(float,[3,size,size,nzer(sum)]);

  prepzernike,size,pupd,cent,cent;    // size = taill tot support, cent=size/2+0.5

  def(,,1) = zernike(4);
  def(,,2) = zernike(2);
  def(,,3) = zernike(3);
  def(,,4:) = prepmirror(size,pupd,nactu,0.09,dmpitch,dmpitch,x0,y0ipupil,ipupil,def(,,1:3),gradG,orthonorm=1) ;

  //building waffle mode
  tmpx= indgen(nactu);
  tmpx = tmpx(,-:1:nactu);
  tmpy = transpose(tmpx);
  waffle = tmpx*0-1;
  x = span(-1,1,nactu)(,-:1:nactu);
  y = transpose(x);
  r = sqrt( x(i,j)^2 + y(i,j)^2 );
  valid = where(r<1.1);
  waffle = waffle^(x+y);

  mywaffle = waffle(valid)(+) * def(,,4:)(,,+);

  nmodes = dimsof(def)(4);
  
  coeff = [gaussdev(nmodes)*.5/(indgen(nmodes))];
  coeff = coeff(*);

  myphase = mywaffle + coeff(+)*def(,,+);

}

/*
//create pupil from adonf imat

tmp = fits_read("~/Downloads/ADONF.fits");
pup = tmp(,,max);
pup1 = float(pup>3500);

require,"/home/brujo/yorick/Ruby/ruby.i";

pup2 = array(0.0,256,256,4);
pup2(,,1) = pup1(1:256,1:256);
pup2(,,2) = pup1(1:256,257:);
pup2(,,3) = pup1(257:,1:256);
pup2(,,4) = pup1(257:,257:);


for (i=2;i<=4;i++) {
res = ruby(crit,image=pup2(,,i),reference=pup2(,,1),variance=1.);
pup2(,,i) = roll(pup2(,,i),[-res(1),-res(2)]);
}

mypup = float(pup2(,,sum)>1.);

//get the center of gravity and center on it
x= float(indgen(256)-256/2-1);
x=x(,-:1:256);
y = transpose(x);
cogx = sum(x*mypup)/sum(mypup);
cogy = sum(y*mypup)/sum(mypup);

mypup2 = roll(mypup,[-cogx,-cogy]);

mypup3=glitch(mypup2,0.7);

thepupil = mypup3>0;


func make_pupil2(dim,pupd,xc=,yc=,real=,cobs=,xcobs=,ycobs=)

{
  if (real == 1) {
    pup = exp(-(mydist(dim,xc=xc,yc=yc)/(pupd/2.))^60.)^0.69314;
  } else {
    //    xc;yc;info,xc;
    //    tv,dist(dim,xc=xc,yc=yc);pause,2000;
    pup = mydist(dim,xc=xc,yc=yc) < (pupd+1.)/2.;
  }
  if (xcobs != []) xc+=xcobs;
  if (ycobs != []) yc+=ycobs;
  if (cobs!=[]) {
    if (real == 1) {
      pup -= exp(-(mydist(dim,xc=xc,yc=yc)/(pupd*cobs/2.))^60.)^0.69314;
    } else {
      pup -= mydist(dim,xc=xc,yc=yc) < (pupd*cobs+1.)/2.;
    }
  }
    
  return pup;
}

pupfit = make_pupil2(256,107,xc=128+0.5,yc=128+0.5,cobs=0.255,xcobs=.5,ycobs=.5);

//so with pixels of 23mas and size = 128, we have pupd = 36
//i.e. a sampling factor of 3.55556
//on this image, we measure a pupd of 107 this leads to a support size, size = 380
//so we put this pupil in a larger support :
thebigpupil = array(0,380,380);
xs = (380-256)/2+1;
ys = (380-256)/2+1;

thebigpupil(xs:xs+256-1,ys:ys+256-1) = thepupil;


 */

/*
stars to remove :
5 6 12 16 21 

valid actus :

ind0=[10,11,12,13,14,15,16,20,21,22,23,24,25,26,27,28,29,30,33,34,35,36,37,38,39,40,41,42,43,44,45,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,139,140,141,142,143,144,145,146,148,149,150,151,152,153,154,155,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,249,250,251,252,253,254,255,256,257,258,259,260,261,264,265,266,267,268,269,270,271,272,273,274,278,279,280,281,282,283,284]

ind4=[15,16,17,18,19,20,27,28,29,30,31,32,33,34,35,36,37,38,43,44,45,46,47,48,49,50,51,52,53,54,55,56,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,361,362,363,364,365,366,367,368,369,370,371,372,373,374,379,380,381,382,383,384,385,386,387,388,389,390,397,398,399,400,401,402]

ind9=[20,21,22,23,24,25,32,33,34,35,36,37,38,39,45,46,47,48,49,50,51,52,53,54,59,60,61,62,63,64,65,66,67,68,69,70,75,76,77,78,79,80,81,82,83,84,85,86,91,92,93,94,95,96,97,98,99,100,101,102,107,108,109,110,111,112,113,114,115,116,117,118,123,124,125,126,127,128,129,130,131,132,133,134,139,140,141,142,143,144,145,146,147,148,149,150,155,156,157,158,159,160,161,162,163,164,170,171,172,173,174,175,176,177,184,185,186,187,188,189]


coeffs=pray_param(41:);
dm0 = coeffs(1:293);
dm4 = coeffs(294:293+416);
dm9 = coeffs(294+416:);

dm0_valid = dm0(ind0);
dm4_valid = dm4(ind4);
dm9_valid = dm9(ind9);

dms_coeffs = _(dm0_valid,dm4_valid,dm9_valid);
  
 */
