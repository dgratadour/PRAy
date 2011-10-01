plug_in,"pray_utils";

require,"string.i";

/*
////////////////////////////////////////////////////////////
////// FROM YAO
////////////////////////////////////////////////////////////

the following routines have been copy-paste from yao
*/

extern _get2dPhase
/* PROTOTYPE
   int _get2dPhase(pointer pscreens, int psnx, int psny, int nscreens, pointer outphase, int phnx, int phny, pointer ishifts, pointer xshifts, pointer jshifts, pointer yshifts)
*/

extern _dist
/* PROTOTYPE
   int _dist(pointer dptr, long dimx, long dimy, float xc, float yc)
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


////////////////////////////////////////////////////////////
////// GENRAL UTILITIES
////////////////////////////////////////////////////////////


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
  y = span(1,-1,nactu)(,-:1:nactu);
  x = span(-1,1,nactu)(-:1:nactu,);
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
        xc = x0 + x(i,j)*dmpitchX;
        yc = y0 + y(i,j)*dmpitchY;
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

func prepcanamir(size,patchDiam,nactu,couplingFact,dmpitchX,dmpitchY,x0,y0,defpup,deftt,norm=)
/* DOCUMENT 
 * prepcanamir(size,patchDiam,nactu,couplingFact,dmpitchX,dmpitchY,x0,y0,defpup,deftt)
 * This routine returns the influence functions in the grid model
 *
 */
{
  if (norm == []) norm = 0;

  nactu = long(nactu(1));
  pup = defpup;
  y = span(1,-1,nactu)(,-:1:nactu);
  x = span(-1,1,nactu)(-:1:nactu,);
  k=0;

  x0       *= patchDiam;
  y0       *= patchDiam;
  x0       += (size/2+0.5);
  y0       += (size/2+0.5);
  dmpitchX *= patchDiam;
  dmpitchY *= patchDiam;

  def = array(0.0, size, size, nactu*nactu);
  gradG = array(0.0, size, size, nactu*nactu,4);
  Ax = -log(couplingFact)/((2./(nactu-1))^2);
  Ay = -log(couplingFact)/((2./(nactu-1))^2);
  
  for(i=1;i<=nactu;i++) {
    for(j=1;j<=nactu;j++) {
      r = sqrt(x(i,j)^2 + y(i,j)^2 );
      if( r<1.1 ) {
        k++;
        xc = x0 - x(i,j)*dmpitchX;
        yc = y0 - y(i,j)*dmpitchY;
        def(,,k) = mygauss2(size,xc,yc,Ax/2.82843,Ay/2.82843,1.,0.,0.);
      }
    }
  }
  def = pup * def(,,:k);

  if (norm==1) {
    nmodes = dimsof(def)(4);
    norm = sum(deftt(,,2)^2);
    
    for (cc=1;cc<=nmodes;cc++) {
      norm_tmp = (sqrt(sum(def(,,cc)^2)/norm))
        def(,,cc) /= norm_tmp; 
    }
  }
  
  return def;
}

func prepbasecano(pup,tt,norm=)
{

  if (norm == []) norm = 0;
  
  valids = where(pup);
  npts = numberof(valids);
  def = pup(,,-::npts-1)*0.;
  for (i=1;i<=npts;i++) def(*,i)(valids(i)) = 1;

  if (norm==1) {
    nmodes = dimsof(def)(4);
    norm = sum(tt(,,2)^2);
    
    for (cc=1;cc<=npts;cc++) {
      norm_tmp = (sqrt(sum(def(,,cc)^2)/norm))
        def(,,cc) /= norm_tmp; 
    }
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

func create_interp_mat(dimx,dimy)
{
  nn = max([dimx,dimy]);
  tmp = indices(nn)(1:dimx,1:dimy,);
  tmp(,,1) -= (dimx/2+1);
  tmp(,,2) -= (dimy/2+1);
  
  mymat = [tmp(,,1)(*)^2,tmp(,,2)(*)^2,tmp(,,1)(*)*tmp(,,2)(*),tmp(,,1)(*),
           tmp(,,2)(*),tmp(,,1)(*)*0+1];
  return transpose(LUsolve(mymat(+,)*mymat(+,))(+,)*mymat(,+));
  // the transpose is to comply with original convention in fitmax
}

func fitmax_multi(im,xmin,xmax,ymin,ymax,interp_mat=,full=)
//same as fitmax but on an arbitrary number of points in x and y
{
  z = im(xmin:xmax,ymin:ymax)(*);
  nx = xmax-xmin+1;
  ny = ymax-ymin+1;

  xm = (xmin+xmax);
  ym = (ymin+ymax);
  if (xm % 2 == 0) {xm /= 2;
  } else xm = xm/2+1;
  if (ym % 2 == 0) {ym /= 2;
  } else ym = ym/2+1;

  if (interp_mat == []) {interp_mat = create_interp_mat(nx,ny);
  } else {if ((dimsof(interp_mat)(sum) != ([2,nx*ny,6])(sum))) {error,"interp_mat has wrong dimensions";}}
    
  p = interp_mat(+,) * z(+); // coeffs de p(1).x^2+p(2).y^2+p(3).xy+p(4).x+p(5).y+p(6)
  A = p(1);
  B = p(2);
  C = p(3);
  denom = C^2.-4*A*B;
  if( denom==0 ) {
    x0 = y0 = 0.0;
  } else {
    x0 = (2.*B*p(4)-p(5)*C)/denom;
    y0 = (2.*A*p(5)-p(4)*C)/denom;
  }
  D = p(6)-(x0*y0*C+A*x0^2.+B*y0^2.);
  x0 += xm;
  y0 += ym;
  if( full==1 ) {
    valmax = D;   // top of the polynom
    if( (B-A)==0 ) {
      t=pi/4;
      if( C==0 ) t=0.;
    } else {
      t = atan(C/(B-A))/2;
    }
    AA = B*sin(t)^2-C*cos(t)*sin(t)+A*cos(t)^2;
    BB = A*sin(t)^2+C*cos(t)*sin(t)+B*cos(t)^2;
    if (D/AA < 0) fwhmx = 1.66*sqrt( -D/2./AA );
    else fwhmx = 0.;
    if (D/BB < 0)fwhmy = 1.66*sqrt( -D/2./BB );
    else fwhmy = 0.;
    
    return [x0,y0,fwhmx,fwhmy,valmax,t];   // t = angle
  } else {
    return [x0,y0];
  }
}

func pup_model(x,a)
{
  // x is the indep var (size of the image)
  // a is the vect of params (10 params)

  pupd = a(1);
  xcp  = a(2);
  ycp  = a(3);
  cob = a(4);
  
  mypup = make_pupil(x(1),pupd,xc=xcp,yc=ycp,cobs=cob,real=1);
  prepzernike,x(1),pupd,xcp,ycp;
  mypup = (a(5)*zernike(1)+a(6)*zernike(2)+a(7)*zernike(3)+
            a(8)*zernike(4)+a(9)*zernike(5)+a(10)*zernike(6))*mypup;

  return mypup;
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

////////////////////////////////////////////////////////////
////// CONTROL INTERFACE
////////////////////////////////////////////////////////////

func pyk_status_push(msg,id=)
{
  if (id==[]) id=1;
  pyk,swrite(format="pyk_status_push(%d,'%s')",id,msg);
}

func null(void) { return; }



////////////////////////////////////////////////////////////
////// DISPLAY FUNCTIONS
////////////////////////////////////////////////////////////

func pray_histeq_scale(z, top=, cmin=, cmax=)
/* DOCUMENT histeq_scale(z, top=top_value, cmin=cmin, cmax=cmax)
     returns a byte-scaled version of the array Z having the property
     that each byte occurs with equal frequency (Z is histogram
     equalized).  The result bytes range from 0 to TOP_VALUE, which
     defaults to one less than the size of the current palette (or
     255 if no pli, plf, or palette command has yet been issued).

     If non-nil CMIN and/or CMAX is supplied, values of Z beyond these
     cutoffs are not included in the frequency counts.

     Identical to histeq_scale except it uses sedgesort instead of sort.
     faster for arrays for which many elements are repeated (e.g.
     CCD arrays where pixels values are integers.
   SEE ALSO: bytscl, plf, pli
 */
{
  if (is_void(top)) top= bytscl([0.,1.])(2);  /* palette size - 1 */
  top= long(top);
  if (top<0 | top>255) error, "top value out of range 0-255";
  y= z(*);
  if (!is_void(cmin)) y= y(where(y>=cmin));
  if (!is_void(cmax)) y= y(where(y<=cmax));
  y= sedgesort(y);
  x= span(0.,1., numberof(y));
  xp= span(0.,1., top+2);
  bins= interp(y, x, xp);
  list= where(bins(dif)<=0.0);
  if (numberof(list)) {
    /* some value (or values) of z are repeated many times --
       try to handle this by adding a small slope to the sorted y */
    dy= y(0)-y(1);
    if (!dy) dy= 1.0;
    for (eps=1.e-10 ; eps<1000.1 ; eps*=10.) {
      bins= interp(y+eps*dy*x, x, xp);
      list= where(bins(dif)<=0.0);
      if (!numberof(list)) break;
    }
    if (eps>1000.) error, "impossible error??";
  }
  return char(max(min(digitize(z,bins)-2,top),0));
}


func escchar(s)
{
  if (s==[]) return;
  s=streplace(s,strfind("_",s,n=20),"!_");
  s=streplace(s,strfind("^",s,n=20),"!^");
  return s;
}


func pray_pltitle(title,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = pltitle_height;
  pltitle_height = long(pltitle_height*defaultdpi/83.);
  
  port= viewport();
  plth=pltitle_height;

  plt, escchar(title), port(zcen:1:2)(1), port(4)+0.005,
    font=pltitle_font, justify="CB", height=plth;

  pltitle_height = plth_save;
}


func pray_xytitles(xtitle,ytitle,adjust,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = pltitle_height;
  pltitle_height = long(pltitle_height*defaultdpi/83.);
  
  curw=current_window();
  if (adjust==[]) {
    adjust = [0.012,0.019];
  }
  xytitles,escchar(xtitle),escchar(ytitle),adjust;

  pltitle_height = plth_save;
}

/*
copy to styc

cp ~/yorick/PRAy/trunk/yorick/Makefile /home/brujo/yorick/canary/styc-nono/yorick/.
cp ~/yorick/PRAy/trunk/yorick/pray_core.i /home/brujo/yorick/canary/styc-nono/yorick/.
cp ~/yorick/PRAy/trunk/yorick/pray.i /home/brujo/yorick/canary/styc-nono/yorick/.
cp ~/yorick/PRAy/trunk/yorick/pray_utils.c /home/brujo/yorick/canary/styc-nono/yorick/.
cp ~/yorick/PRAy/trunk/yorick/pray_utils.i /home/brujo/yorick/canary/styc-nono/yorick/.
cp ~/yorick/PRAy/trunk/yorick/widget_pray.i /home/brujo/yorick/canary/styc-nono/yorick/.

cp ~/yorick/PRAy/trunk/widgets/widget_pray.py /home/brujo/yorick/canary/styc-nono/widgets/.
cp ~/yorick/PRAy/trunk/glade/widget_pray.glade /home/brujo/yorick/canary/styc-nono/glade/.


 */
