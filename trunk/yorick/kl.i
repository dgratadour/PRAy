/*
 Collection of routines for building the KL of a Kolmogorov statistics
 Derived from a collection of routines developped both at ESO and ONERA
 The main routine is the last one ...

 for a Kolmogorov statistics :
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64);
 or
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="kolmo");


 for a Von Karman statistics :
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman");
 or
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman",outscl=3);
 default is : an outter scale of 3 times the size of the telescope

 res is a 128x128x150 array containing the 150 first KL
 of a kolmogorov or Von Karman stat

 D. Gratadour Feb 2006
*/

require,"digit2.i";
require,"pray_core.i";

func dblindgen(n)
 /*DOCUMENT res=dblindgen(size)

 D. Gratadour Feb 2006

 This routine returns a size x size array containing increasing indices
 from 1 to size x size.

 SEE ALSO : indgen, indices
  */
{
 n=long(n);
 return reform(indgen(n*n),[2,n,n]);
}

func polar_coord(r,&mask,&rho,&phi,&pts,occ=,xcent=,ycent=,\
                verbose=,leq=,btw4pix=,sizemin=,dbprec=)
 /* DOCUMENT polar_coord,radius,mask,rho,phi;
          or polar_coord,radius,mask,rho,phi,pts,sizemin=1;

 D. Gratadour Jan 2006

 Calculation of polar coordinates rho and phi and an intensity mask
 (the pupil) of a telescope of radius r
 Derived from an IDL routine (polaire2.pro) written by L. Mugnier

 INPUTS :
 r = the radius of the mask

 OUTPUT :
 mask = the intensity mask (2d image)
 rho  = the width coordinate (2d image)
 phi  = the angle coordinate (rad) (2d image)
 pts  = (optional) indices of valid (non null) points of the
        pupil (1d vector) if flag sizemin is set 

 OPTIONAL :
 occ     = the occultation level (<=1)
 xcent   = the x position of the center point of the mask 
 xcent   = the y position of the center point of the mask 
 verbose = flag to get info on the process (0/1)
 leq     = flag to set the limits of the mask (<= radius or < radius)
 btw4pix = flag to set the center of the mask on 1 pixel or
           in between 4 pixels (0/1)
 sizemin = flag to set the outputs in a minimum size arrays (null points
           are eliminated) (0/1)

 SEE ALSO : ...
 */
{
 if (r==[]) r=128;
 if (!is_set(occ)) occ=0.0;
 if (!btw4pix) diam=long(2*floor(r)+1);
 else diam=long(2*round(r));

 if (!is_set(xcent)) xcent=float((diam-1)/2.);  
 if (!is_set(ycent)) ycent=float((diam-1)/2.);

 if (verbose) {
   write,format="Pupil diameter : %d\n",diam;
   if (xycent!=[]) {
     write,format="Pupil X center : %f\n",xycent(1);
     write,format="Pupil Y center : %f\n",xycent(2);
   }
 }

 x=float(((indices(diam))(,,1)-1.0) % diam);
 y=transpose(x);
 x-=xcent;
 y-=ycent;

 if (dbprec) {
   rho=double(sqrt(x^2+y^2))/double(r);
   phi=double(atan(y,x+(rho==0)));
 }
 else {
   rho=float(sqrt(x^2+y^2))/float(r);
   phi=float(atan(y,x+(rho==0)));
 }

 if (leq) mask=((rho<=1) & (rho>=occ));
 else {
   if (dbprec) mask=double((rho<1) & (rho>=occ));
   else mask=float((rho<1) & (rho>=occ));
 }

 if (sizemin) {
   if (verbose) write,"We will keep only the non-null points";
   pts=where(mask != 0.);
   if (dimsof(pts) != []) {
     sizepts=(dimsof(pts))(2);
     if (verbose) write,"Number of null points : ",sizepts;
     rho=rho(pts);
     phi=phi(pts);
     mask=mask(pts);
   }
   else {
     if (verbose) write,"All points are non-null !";
   }
 }
}


struct gkl_basis_struct
{
 long   nr;            
 long   ni;  
 long   np;       
 long   nfunc;   
 float   ri;  
 pointer   radp;   
 pointer evals;   
 long   nord;   
 pointer   npo;   
 pointer   ord;   
 pointer   rabas;   
 pointer   azbas;   
};

struct geom_struct
{
 pointer   px;            
 pointer   py;  
 pointer   cr;       
 pointer   cp;   
 pointer   pincx;  
 pointer  pincy;   
 pointer pincw;   
 pointer   ap;   
 long   ncp;   
 long   ncmar;   
};


func radii(nr,np,ri)
 /*DOCUMENT res=radii(NumberOfR,NumberOfPhi,Dim)

 D. Gratadour Feb 2006

 This routine generates an nr x np array with np copies of the
 radial coordinate array. Radial coordinate span the range from
 r=ri to r=1 with successive annuli having equal areas (ie, the
 area between ri and 1 is divided into nr equal rings, and the
 points are positioned at the half-area mark on each ring). There
 are no points on the border.     

 SEE ALSO : polang
  */
{
 r2 = ri^2 +(float(indgen(nr)-1)+0.)/nr*(1.0 - ri^2);
 rs = sqrt(r2);
 return rs*array(1.,np)(-,);
}

func polang(r)
 /*DOCUMENT res=polang(RadialCoordArray)

 D. Gratadour Feb 2006

 This routine generates an array with the same dimensions as r,
 but containing the azimuthal values for a polar coordinate system.     

 SEE ALSO : radii
  */
{
 s =  dimsof(r);
 nr = s(2);
 np = s(3);
 phi1 = float(indgen(np)-1)/float(np)*2.*pi;
 return phi1(-,)*array(1.,nr);
}

func setpincs(ax,ay,px,py,ri,&pincx,&pincy,&pincw)
 /*DOCUMENT res=polang(RadialCoordArray)

 D. Gratadour Feb 2006

 This routine determines a set of squares for interpolating
 from cartesian to polar coordinates, using only those points
 with ri < r < 1     

 SEE ALSO : pcgeom
  */
{
 s = dimsof(ax);
 nc = s(2);
 s = dimsof(px);
 nr = s(2);
 np = s(3);
 dcar = (ax(nc) - ax(1)) / (nc-1);
 ofcar = ax(1,1);

 rlx = (px - ofcar)/dcar;
 rly = (py - ofcar)/dcar;
 lx = long(rlx);
 ly = long(rly);
 shx = rlx - lx;
 shy = rly - ly;

 pincx=[lx,lx+1,lx+1,lx]+1;
 pincy=[ly,ly,ly+1,ly+1]+1;

 pincw=[(1-shx)*(1-shy),shx*(1-shy),shx*shy,(1-shx)*shy];

 axy = ax^2 + ay^2;
 axyinap = clip(axy,ri^2.+1.e-3,0.999);
 sizeaxyinap=(dimsof(axyinap))(2);
 pincw = pincw*axyinap(pincx+(pincy-1)*sizeaxyinap);
 pincw = pincw*(1.0/pincw(,,sum))(,,-);
}

func pcgeom (nr,np,ncp,ri,ncmar,ap)    
 /*DOCUMENT geom=pcgeom(nr, np, ncp, ri, ncmar,ap)

 D. Gratadour Feb 2006

 This routine builds a geom_struct. px and py are the x and y
 coordinates of points in the polar arrays.  cr and cp are the
 r and phi coordinates of points in the cartesian grids. ncmar
 allows the possibility that there is a margin of ncmar points
 in the cartesian arays outside the region of interest


 SEE ALSO : setpincs, set_pctr
  */
{
 nused = ncp - 2*ncmar;
 ff = 0.5 * nused;
 hw =  float(ncp-1)/2;

 r = radii(nr,np,ri); 
 p = polang(r);

 px0 = r * cos(p);
 py0 = r * sin(p);
 px = ff * px0 + hw;
 py = ff * py0 + hw;
 ax = float(dblindgen(ncp)-1) % ncp - 0.5 * (ncp-1);
 ax = ax / (0.5 * nused); 
 ay = transpose(ax);

 setpincs, ax, ay, px0, py0, ri,pincx, pincy, pincw;
 dpi = 2 * pi;
 cr2 = (ax^2 + ay^2); 
 ap = clip(cr2,ri^2+1.e-3,0.999);
 //cr = (cr2 - ri^2) / (1 - ri^2) * nr - 0.5; 
 cr = (cr2 - ri^2) / (1 - ri^2) * nr; 
 cp = (atan(ay, ax) + dpi) % dpi;
 cp = (np / dpi) * cp;

 cr = clip(cr,1.e-3,nr-1.001);
 //fudge -----, but one of the less bad ones
 cp = clip(cp,1.e-3,np -1.001);
 //fudge -----  this is the line which
 //gives that step in the cartesian grid
 //at phi = 0.

 geom = geom_struct();
 geom.px=&px;
 geom.py=&py; 
 geom.cr=&cr;
 geom.cp=&cp;
 geom.pincx=&pincx;
 geom.pincy=&pincy;
 geom.pincw=&pincw;
 geom.ap=&ap;
 geom.ncp=ncp;
 geom.ncmar=ncmar; 

 return geom;
}

func set_pctr(bas, ncp =, ncmar=)
 /*DOCUMENT geom=set_pctr(bas, ncp =, ncmar=)

 D. Gratadour Feb 2006

 This routine calls pcgeom to build a geom_struct with the
 right initializations. bas is a gkl_basis_struct built with
 the gkl_bas routine.

 SEE ALSO : pcgeom, setpincs, gkl_bas
  */
{
 if (!is_set(ncmar)) ncmar = 2;
 if (!is_set(ncp)) ncp = 128;

 return pcgeom(bas.nr,bas.np,ncp,bas.ri,ncmar,ap);
}

func pol2car(cpgeom,pol,mask=)
 /*DOCUMENT cart=pol2car(cpgeom, pol, mask=)

 D. Gratadour Feb 2006

 This routine is used for polar to cartesian conversion.
 pol is built with gkl_bas and cpgeom with pcgeom.
 However, points not in the aperture are actually treated
 as though they were at the first or last radial polar value
 -- a small fudge, but not serious  ?*******

 SEE ALSO : pcgeom, gkl_bas
  */
{
 if (sae) error;
 cd = bilinear(pol, *cpgeom.cr+1, *cpgeom.cp+1);
 if (mask!=[]) cd = cd*(*cpgeom.ap);
 return cd;
} 

func kolstf(dvec)
 /*DOCUMENT var=kolstf(dvec)

 D. Gratadour Feb 2006

 This routine returns the kolmogorov phase variance at spatial
 dimension (inverse of the spatial frequency) dvec

 SEE ALSO : 
  */
{
 return  6.88 * dvec^(5./3.);
}

func karmanstf(dvec,outscl=)
 /*DOCUMENT var=kolstf(dvec)

 D. Gratadour Feb 2006

 This routine returns the Von Karman phase variance at spatial
 dimension (inverse of the spatial frequency) dvec. Same as kolstf
 but with a correcting factor to account for the outter scale.
 The latter should be in units of telescope diameter

 SEE ALSO : 
  */
{
 if (dimsof(outscl)==[]) outscl = 3.;
 return 6.88 * dvec^(5./3.)*(1-1.485*(dvec/outscl)^(1./3.)+\
                             5.383*(dvec/outscl)^(2)-6.281*\
                             (dvec/outscl)^(7./3.));
}

func gkl_radii(nr,ri)
 /*DOCUMENT rad=gkl_radii(nr,ri)

 D. Gratadour Feb 2006

 This routine generates an array of radial polar coordinates along
 which the KL are generated. nr is the number of elements and ri is
 the maximum radius.

 SEE ALSO : 
  */
{
 d = (1.-ri*ri)/nr;
 //    rad2 = ri^2 + d/2. + d * float(indgen(nr)-1);
 //    rad2 = ri^2 + d * float(indgen(nr)-1);
 rad2 = ri^2 +d/16.+ d * float(indgen(nr)-1);
 //  rad2 = ri^2 +d/14.+ d * float(indgen(nr)-1);  // nr=64,128
 //  rad2 = ri^2 +d/10.+ d * float(indgen(nr)-1);
 rad = sqrt(rad2);

 return rad;
}

func gkl_mkker(ri,nr,rad,funct=,outscl=)
 /*DOCUMENT 

 D. Gratadour Feb 2006

 This routine generates the kernel used to find the KL modes.
 The  kernel constructed here should be simply a discretization
 of the continuous kernel. It needs rescaling before it is treated
 as a matrix for finding  the eigen-values. The outter scale
 should be in units of the diameter of the telescope.

 SEE ALSO : 
  */
{
 nth = 5*nr;
 kers  = array(float,[3,nr, nr, nth]);
 cth = cos(float(indgen(nth)-1)*(2.*pi/nth));
 dth = 2.*pi/nth;
 fnorm = -1./(2*pi*(1.-ri^2))*0.5;
 //the 0.5 is to give  the r^2 kernel, not
 //the r kernel
 for (i =1;i<=nr;i++) { 
   for (j=1;j<=i;j++) {
     te = 0.5*sqrt(rad(i)^2+rad(j)^2-(2*rad(i)*rad(j))*cth);
     //te in units of the diameter, not the radius
     if (funct=="kolmo") te = kolstf(te);
     if (funct=="karman") te = karmanstf(te,outscl=outscl);
     if ((funct!="kolmo") & (funct!="karman")) {
       write,"The statistics is not known !";
       error;
     }
     kelt =  fnorm * dth * float (fft(te,-1));
     kers (i, j,) = kelt;
     kers (j, i,) = kelt;
   }
   if (is_set(verbose)) write, i;
 }
 if (is_set (verbose))  write," ";

 return kers;

}

func piston_orth(nr)
{
 s = array(float,[2,nr,nr]);
 for (j=1;j<=nr-1;j++) {
   rnm = 1./sqrt (float((j)*(j+1)));
   s(1:j,j) = rnm;
   s(j+1,j)= -1*(j)*rnm;
 }
 rnm = 1./sqrt (nr);
 s(,nr) = rnm;
 return s;
}

func gkl_fcom(kers,ri,nf,&evals,&nord,&npo,&ord,&rabas)
 /*DOCUMENT 

 D. Gratadour Feb 2006

 This routine does the work : finding the eigenvalues and
 corresponding eigenvectors. Sort them and select the right
 one. It returns the KL modes : in polar coordinates : rabas
 as well as the associated variance : evals. It also returns
 a bunch of indices used to recover the modes in cartesian
 coordinates (nord, npo and ord).

 SEE ALSO : gkl_bas
  */
{
 s = dimsof(kers);
 nr = s(2);
 nt = s(4);
 nxt = 1;
 fktom =  (1.-ri^2)/nr;
 fevtos = sqrt(2*nr);
 evs = array(float,[2,nr,nt]);
 //ff isnt used - the normalisation for
 //the eigenvectors is straightforward:
 //integral of surface^2 divided by area = 1,
 //and the cos^2 term gives a factor
 //half, so multiply zero order by
 //sqrt(n) and the rest by sqrt (2n)

 //zero order is a special case...
 //need to deflate to eliminate infinite eigenvalue - actually want
 //evals/evecs of zom - b where b is big and negative
 zom = kers(,,1);
 s = piston_orth(nr);
 ts =transpose(s);
 b1 = ((ts(,+)*zom(+,))(,+)*s(+,))(1:nr-1, 1:nr-1);

 newev = SVdec(fktom*b1,v0,vt);

 v1 = array(float,[2,nr, nr]);
 v1(1:nr-1,1:nr-1) = v0;
 v1(nr,nr) = 1;

 vs = s(,+)*v1(+,);
 grow,newev,0;
 evs(,nxt) = float(newev);
 kers (,, nxt) = sqrt(nr)*vs;
 // the rest are more straightforward
 nxt = 2;
 do {
 newev = SVdec(fktom*kers(,,nxt),vs,vt);
 evs(,nxt) = float(newev);
 kers (,,nxt) = sqrt(2.*nr)*vs;
 mxn = max(float(newev));
 egtmxn = floor(evs(, 1:nxt)>mxn);
 nxt = nxt + 1;
 } while ((2*sum(egtmxn)-sum(egtmxn(,1))) < nf);
 nus = nxt - 1;

 kers = kers (,,1:nus);
 evs = reform (evs (, 1:nus), nr*nus);
 a = (sort(-1.*evs))(1:nf);
 //every eigenvalue occurs twice except
 //those for the zeroth order mode. This
 //could be done without the loops, but
 //it isn't the stricking point anyway...
 no = 1;
 ni = 1;
 oind = array(long,nf+1);
 do {
      if (a(ni) < nr+1) {
        oind(no) = a(ni);
        no = no + 1;
      } else {
        oind(no) = a(ni);
        oind(no+1) = a(ni);
        no = no + 2;
      }
      ni = ni + 1;
 } while (no < (nf+1));

 oind = oind (1:nf);
 tord = (oind-1)/nr+1;
 odd = ((long(indgen(nf)-1) % 2) == 1);
 pio = (oind-1) % nr +1;

 evals = evs(oind);
 ord = 2 *(tord-1) - floor(tord>1 & (odd))+1;

 nord = max(ord);
 rabas = array(float,[2,nr, nf]);
 sizenpo=long(max(ord));
 npo = array(long,sizenpo);

 for (i=1;i<=nf;i++) {
   npo(long(ord(i))) = npo(long(ord(i))) + 1;
   rabas(, i) = kers (, pio(i), tord(i));
 }
}

func gkl_mkazi(nord, np)
{
 gklazi = array(float,[2,long(1+nord), np]);
 th = float(indgen(np)-1)*(2.*pi/ np);

 gklazi (1,) = 1.0;
 for (i = 2; i<=nord;i+=2)  gklazi (i,) = cos (((i-1)/2+1) * th);
 for (i = 3; i<=nord;i+=2)  gklazi (i,) = sin (((i-1)/2) * th);
 return gklazi;
}

func gkl_bas(ri=,nr=,np=,nfunc=,verbose=,funct=,outscl=)
 /*DOCUMENT 

 D. Gratadour Feb 2006

 This routine uses the output of gkl_fcom to fill the gkl_base_struct.

 SEE ALSO : gkl_fcom
  */
{
 if (!is_set(ri)) ri = 0;
 if (!is_set(nr)) nr = 40;
 if (!is_set(np)) np = long(5*nr);
 if (!is_set(nfunc)) nfunc = 500L;

 nr = long(nr);
 np = long(np);

 if ((nr * np)/ nfunc < 8) {
   if (is_set(verbose)) write,"warning: you may need a finer ",\
                          "radial sampling ";
   if (is_set(verbose)) write, "(ie, increased nr) to generate ",\
                          nfunc, "  functions";
 } else if ((nr * np)/ nfunc > 40) {
   if (is_set(verbose)) write,"note, for this size basis ",\
                          "radial discretization on ", nr;
   if (is_set(verbose)) write, "points is finer than necessary",\
                          "-it should work, but you ";
   if (is_set(verbose)) write, "could take a smaller nr without",\
                          "loss of accuracy";
 }


 radp = gkl_radii(nr, ri);

 kers = gkl_mkker(ri, nr, radp,funct=funct,outscl=outscl);

 gkl_fcom,kers,ri,nfunc,evals,nord,npo,ord,rabas;

 azbas = gkl_mkazi(nord, np);

 gklbasis = gkl_basis_struct();
 gklbasis.nr=nr;
 gklbasis.np=np;
 gklbasis.nfunc=nfunc; 
 gklbasis.ri=ri;
 gklbasis.radp=&radp;
 gklbasis.evals=&evals;
 gklbasis.nord=nord;
 gklbasis.npo=&npo;
 gklbasis.ord=&ord;
 gklbasis.rabas=&rabas;
 gklbasis.azbas=&azbas;

 return gklbasis;
}

func gkl_sfi(bas, i)
 /*DOCUMENT 

 D. Gratadour Feb 2006

 This routine returns the i'th function from the generalised KL
 basis bas. bas must be generated first with gkl_bas.

 SEE ALSO : gkl_bas
  */
{    
 if (i>bas.nfunc) { 
   write, "the basis only contains ", nfunc, "functions";
   return 0;
 }

 nr = bas.nr;
 np = bas.np;
 ordp = *bas.ord;
 ord=long(ordp(i));

 rabasp=*bas.rabas;
 rabas=rabasp(,i);

 azbasp=*bas.azbas;
 azbas=azbasp(ord, );

 sf1=array(double,[2,nr,np]);
 sf1(,*)=rabas;

 sf2=array(float,[2,np,nr]);
 sf2(,*)=azbas;  

 sf = sf1*transpose(sf2);
 return sf;
}


func make_kl(nmax,dim,&var,&outpolarbase,&pupil,oc=,nr=,nopup=,\
            funct=,outscl=,verbose=)
/* DOCUMENT 
 for a Kolmogorov statistics :
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64);
 or
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="kolmo");


 for a Von Karman statistics :
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman");
 or
 res=make_kl(150,128,varkl,outbas,pup1,oc=0.12,nr=64,funct="karman",outscl=5);

 the outter scale is in units of the telescope diameter
 default is : an outter scale of 3 times the size of the telescope

 D. Gratadour Feb 2006

 This routine is the main program. It returns nmax generalized
 KL in an array dim x dim x nmax. It also returns the associated
 variance as well as the pupil and the polar base used for their
 calculation. Optional keywords includes any occultation, the
 number of samples for the radial coordinate and a flag to avoid
 pupil multiplication.

 SEE ALSO : polar_coord, gkl_bas, set_pctr
*/
{
 if (pupil==[]) polar_coord,dim/2.,pup,rho,phi,occ=oc,btw4pix=1;
 else pup=pupil;

 if (!is_set(nr)) nr=64;

 if (dimsof(funct)==[]) {
   write,"using the Kolmogorov model";
   funct="kolmo";
 }

 polarbase = gkl_bas(ri=oc,nr=nr,np=(2*pi*nr),nfunc=nmax,\
                     funct=funct,outscl=outscl,verbose=verbose);

 outpolarbase = polarbase;

 pc1 = set_pctr(polarbase, ncp= dim);

 kl = array(float,[3,long(dim),long(dim),nmax]);

 if (is_set(nopup)) {
   for (i=1;i<=nmax;i++) kl(,,i)=pol2car(pc1, gkl_sfi(polarbase,i));
 } else {
   for (i=1;i<=nmax;i++) {
     sae=0;
     //      if (i==8) sae=1;
     kl(,,i)=pol2car(pc1, gkl_sfi(polarbase,i))*pup;
   }
 }

 pupil =  pup; 
 var =  *polarbase.evals;

 return kl;
}


