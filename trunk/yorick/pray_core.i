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
  pointer pupmap;
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
  long    filt_quad;
};


/*
 ____  _                      ____  _                    _ _         
|  _ \| |__   __ _ ___  ___  |  _ \(_)_   _____ _ __ ___(_) |_ _   _ 
| |_) | '_ \ / _` / __|/ _ \ | | | | \ \ / / _ \ '__/ __| | __| | | |
|  __/| | | | (_| \__ \  __/ | |_| | |\ V /  __/ |  \__ \ | |_| |_| |
|_|   |_| |_|\__,_|___/\___| |____/|_| \_/ \___|_|  |___/_|\__|\__, |
                                                               |___/ 
*/

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

func pray_coeff2psfs(coeff,&ampliPup,&ampliFoc,poffset=,pupmap=)
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
    if (pupmap != [])
      msk = ipupil * pupmap;
    else msk = ipupil;
    psfs(,,i) = phase2psf(tmp(_n1:_n2,_n1:_n2),msk(_n1:_n2,_n1:_n2),\
                          size,amp1,amp2,fftws=fftws);
    ampliPup(,,i) = amp1;
    ampliFoc(,,i) = amp2;
  }
  return psfs;
}

func pray_init(xpos,ypos,tmodes=,mir_params=,tiptilt=,poffset=)
/* DOCUMENT pray_init
 *  
 * pray_init,xpos,ypos,tmodes=,mir_params=,tiptilt=,poffset=
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
        /*
        if (filter_quad == 0) {
          for (i=1;i<=nzer(k);i++) {
            cpt ++;
            def(,,cpt) = zernike(i+1);
            //def(,,cpt) = zernike(selz(i));
          }
        } else {
        */
          for (i=7;i<=nzer(k)+6;i++) {
            cpt ++;
            def(,,cpt) = zernike(i);
          }
          //}
      }
    }

    if(tmodes == 1) {
      pup1 = [];
      if (k==1) {
        foc = zernike(4);
        kl = make_kl(nzer(k),pupd,v,obas,ipupil,oc=cobs,nr=128,nopup=1);
        //    kl = order_kls(kl,patchDiam,upto=20);
        def(size/2-pupd/2+1:size/2+pupd/2,            \
            size/2-pupd/2+1:size/2+pupd/2,2:nzer(k)) = kl(,,:-1);
        def *= ipupil;

        // filtering focus : not working
        //res=def-foc(,,-::nzer(k)-1)(,,+)*(diag((def(*,)(where(ipupil),)(+,)*foc(*)(where(ipupil))(+))))(+,);

        tmp = def(,,2:3);
        def(,,1:2) = tmp;
        def(,,3) = foc;
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
      // best fit :
      if (mir_params == []) {
        pxm = 0.485998;
        pym = 0.499824;
        x0m = 0.00156321;
        y0m = 0.00194297;
      } else {
        pxm = mir_params(1);
        pym = mir_params(2);
        x0m = mir_params(3);
        y0m = mir_params(4);
      }
      //deftt = [zernike(2),zernike(3),zernike(4)];
      tmp = prepcanamir(size, pupd,ceil(sqrt(nzer(1))),0.2,pxm,pym,x0m,y0m,ipupil);
      def(,,1:3) = [zernike(2),zernike(3),zernike(4)];
      def(,,4:nzer(1)) = tmp;
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
  pupmap      = *extra.pupmap;
  if (numberof(pupmap)==1) pupmap = [];
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
    def       = *(pray_data.def);
    x0        = param(1);
    y0        = param(2);
    dmpitchX  = param(3);
    dmpitchY  = param(4);
    gradG     = 0.0;
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
    norm_im = images(,,avg,1)(*)(sum);
    new_obj /= norm_obj;
    new_obj *= norm_im;
    grad_obj /= norm_obj;
    grad_obj *= norm_im;
    ftobject = fft(new_obj,1);
  }
  
  //-----------------------------------------------------------------
  // First dealing with in-focus images
  tmp = coeffs;
  psfs = pray_coeff2psfs(tmp,ampliPup,ampliFoc,poffset=poffset,pupmap=pupmap);
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
       
    psfs = pray_coeff2psfs(tmp,ampliPup,ampliFoc,poffset=poffset,pupmap=pupmap);
 
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
 

func pray(images,xpos,ypos,deltaFoc,sigma,object,lambda,nzer=,alt=,teldiam=,cobs=,osampl=,disp=,verbose=,threshold=,nbiter=,tmodes=,tiptilt=,variance=,guess=,scale=,diff_tt=,fit_starpos=,fit_object=,pup_params=,mir_params=,script=)
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
  if (!is_set(disp)) disp=0;
  if (!is_set(script)) script = 0;
                               
  // sizes init and check
  dims = dimsof(images);
  if (dims(1) != 4) {
    if (pray_gui) 
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","images format must follow : [size,size,ntarget,ndefoc]");
    else
      write,"images format must follow : [size,size,ntarget,ndefoc]"
    return 1;
  }
  size = dims(2);
  ntarget = dims(4);
  ndefoc = dims(5)-1;
  if ((numberof(xpos) != ntarget) || (numberof(ypos) != ntarget)) {
    if (pray_gui)
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","incompatible images and xpos size");
    else
      write,"incompatible images and xpos size"
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

  // checking if there are a pupil map parameters
  if (pup_params != []) {
    if (numberof(pup_params) != 7) error,"wrong pup_params size";
    a = _(float(pupd),size/2+0.5,size/2+0.5,pup_params);
    tmp = pup_model(size,a);
    pray_data.pupmap = &float(tmp/max(tmp));
    cobs = pup_params(1);
  } else pray_data.pupmap = &([0.0f]);

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
  if (pray_gui)
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Initializing the PRAy workspace ...");
  else
    write,"Initializing the PRAy workspace ..."
  
      pray_init,xpos,ypos,tmodes=tmodes,mir_params=mir_params;

  if (pray_gui)
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Done ...");
  else
    write,"Done ...";

  // ... testing the size and content of variance
  if (variance==[]) {
    variance = sigma^2;
    sz_variance = 1;
  } else {
    sz_variance = numberof(variance);
    if ((sz_variance != 1) & (sz_variance != numberof(deltaFoc)) & (sz_variance != numberof(images))) {
      if (pray_gui)
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "variance must be a scalar or of same size than image.");
      else
        write,"variance must be a scalar or of same size than image.";
      return 1;
    } else {
      if (min(variance) <= 0.) {
        if (pray_gui)
          pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')", \
                          "variance must be strictly positive");
        else
          write,"variance must be strictly positive";
        return 1;
      }
    }
  }
  
  // ... testing guess and create one if nill
  norm = 0.;//guess*0.;
  if (guess != []) {
    if (numberof(guess) != nzer(*)(sum)) {
      if (pray_gui)
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "Guess does not have the proper size. Initializing to 0s.");
      else
        write,"Guess does not have the proper size. Initializing to 0s.";
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
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",\
               cc,pname);
  }
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_active(%d)",0);
  pray_selected_error=1;
  */
  if (pray_gui) {
    pray_update_zernike_table,pray_param,pray_data.lambda,nzer(1);
    pyk,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_active(%d)",0);
    pray_selected_display=1;
  }

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

  if (!script) pray_loop;
  else pray_loop_script;
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
      if (pray_gui) {
        pray_progressbar_text,"I reached the max number of iterations";
        pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                        "End of The Minimization process");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "---------------------------------------------------------------");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "I reached the max number of iterations");
      } else {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of iterations";
      }
      stop_pray_loop = 1;
    }

    if (pray_eval == pray_data.neval) {
      if (pray_gui) {
        pray_progressbar_text,"I reached the max number of eavaluation";
        pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                        "End of The Minimization process");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "---------------------------------------------------------------");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",   \
                        "I reached the max number of evaluations");
      } else {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the max number of evaluations";
      }
      stop_pray_loop = 1;
    }
    if ((dispok) && (pray_gui)) {
      pray_disp_select;
      pray_disp_error;
    }

    if (pray_gui) {
      pray_progressbar_frac,float(pray_iter)/pray_data.nbiter;
      pray_progressbar_text,swrite(format="iter %d of %d",pray_iter,pray_data.nbiter);
    }
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
    if (pray_gui) {
      pray_progressbar_frac,1.0;
      pray_progressbar_text,"I reached the convergence threshold";
      pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')","End of The Minimization process");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","I reached the convergence threshold");
    } else {
      write,"End of The Minimization process";
      write,"---------------------------------------------------------------";
      write,"I reached the convergence threshold";
    }
    stop_pray_loop = 1;
  }    
 
  if (pray_gui) pray_update_zernike_table,pray_param,pray_data.lambda,(*pray_data.nzer)(1);

  if (stop_pray_loop) {
    stop_pray_loop=0;
    newstring = swrite(format="Value of the criterion after %d iterations : %23.2e",pray_iter,pray_error(pray_param,grad_fin,pray_data));
    if (pray_gui) {
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
      pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
      
      pray_disp_select;
      pray_disp_error;
      pyk,cmd_pray+"y_on_loop_finished()";
    } else {
      write,"---------------------------------------------------------------";
      write,newstring;
    }
    return;
  }
  
  if (!one) {
    after,0.1,pray_loop;
  }
}

func pray_loop_script(void)
{
  extern pray_data,pray_iter;
  extern pray_param,dispok;
  extern stop_pray_loop;
  extern elapsed,cpu_start,pray_task,pray_eval,pray_step,pray_wfs,pray_f,pray_g;

  while (!stop_pray_loop) {
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
        if (pray_gui) {
          pray_progressbar_text,"I reached the max number of iterations";
          pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                     "End of The Minimization process");
          pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",      \
                     "---------------------------------------------------------------");
          pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",      \
                     "I reached the max number of iterations");
        } else {
          write,"End of The Minimization process";
          write,"---------------------------------------------------------------";
          write,"I reached the max number of iterations";
        }
        stop_pray_loop = 1;
      }
      
      if (pray_eval == pray_data.neval) {
        if (pray_gui) {
          pray_progressbar_text,"I reached the max number of eavaluation";
          pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')", \
                     "End of The Minimization process");
          pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",      \
                     "---------------------------------------------------------------");
          pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",      \
                     "I reached the max number of evaluations");
        } else {
          write,"End of The Minimization process";
          write,"---------------------------------------------------------------";
          write,"I reached the max number of evaluations";
        }
        stop_pray_loop = 1;
      }
      if ((dispok) && (pray_gui)) {
        pray_disp_select;
        pray_disp_error;
      }
      
      if (pray_gui) {
        pray_progressbar_frac,float(pray_iter)/pray_data.nbiter;
        pray_progressbar_text,swrite(format="iter %d of %d",pray_iter,pray_data.nbiter);
      }
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
      if (pray_gui) {
        pray_progressbar_frac,1.0;
        pray_progressbar_text,"I reached the convergence threshold";
        pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')","End of The Minimization process");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
        pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","I reached the convergence threshold");
      } else {
        write,"End of The Minimization process";
        write,"---------------------------------------------------------------";
        write,"I reached the convergence threshold";
      }
      stop_pray_loop = 1;
    }    
    
    if (pray_gui) pray_update_zernike_table,pray_param,pray_data.lambda,(*pray_data.nzer)(1);
  }
  
  stop_pray_loop=0;
  newstring = swrite(format="Value of the criterion after %d iterations : %23.2e",pray_iter,pray_error(pray_param,grad_fin,pray_data));
  if (pray_gui) {
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","---------------------------------------------------------------");
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
    
    pray_disp_select;
    pray_disp_error;
    pyk,cmd_pray+"y_on_loop_finished()";
  } else {
    write,"---------------------------------------------------------------";
    write,newstring;
  }
  return;
}

func pray_update_comments(iter, eval, cpu,f,gnorm,step)
{
  extern cmd_pray;

  newstring = swrite(format=" %5d %5d %10.3f     %+-10.6e  %-9.1e  %-9.1e",iter, eval, cpu,f , gnorm, step);
  if (pray_gui) {
    pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')"," ITER    EVAL     CPU [s]          FUNC          GNORM   STEPLEN  ");
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","------    ------     ----------   -----------------------    ---------     ---------");
    pyk,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
  } else {
    if (iter == 1) {
      write," ITER    EVAL     CPU [s]        FUNC       GNORM     STEPLEN  ";
      write,"------  ------  ----------   ----------  -----------    -------";
      write,newstring;
    } else write,format="\r %s \n",newstring;
  }
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
  //pyk_pray,swrite(format=cmd_pray+"glade.get_widget('deltafoc_entry').set_text('%s')",box_content);ic
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

func stop_pray(void)
{
  extern stop_pray_loop;

  stop_pray_loop = 1;
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
