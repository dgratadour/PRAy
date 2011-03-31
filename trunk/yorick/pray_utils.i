/*
Most of these lines of code have been inspired by previous work from francois rigaut
for the gemini mcao system user interface (myst)
 */


require,"string.i";

func pyk_status_push(msg,id=)
{
  if (id==[]) id=1;
  pyk,swrite(format="pyk_status_push(%d,'%s')",id,msg);
}

func null(void) { return; }

func gui_update(void)
{
  //pyk_status_push,"STyC version "+styc_version+" ready !";  
  return;
}

func styc_log(text,close=)
{
  extern stycLogfile,dateTime;

  if (stycLogfile==[]) {
    dateTime = styctop+"/logs/styc_"+styc_get_timestamp()+".log";
    logfile = open(dateTime,"w");
  }
  
  if (close) {
    close,logfile;
    logfile=[];
    return;
  }
  
  if (text==[]) return;

  write,logfile,gettime()+" "+text;
  fflush,logfile;
}

func styc_get_timestamp(&time,&date)
{
  time = gettime();
  date = getdate();
  return sum(strpart(date,[[0,2],[3,5],[6,8]]))+"_"+ \
         sum(strpart(time,[[0,2],[3,5],[6,8]]));
}

func styc_parse_flags(args)
{
  extern styc_simul;

  if (numberof(args)<4) return;
  
  flags = args(4:);
  nflags = numberof(flags);
  valid = array(0,nflags);

  for (i=1;i<=nflags;i++) {
    if (flags(i)=="--simul") {
      styc_simul=1;
      valid(i) = 1;
    }
    if (flags(i)=="--load") {
      if (i==nflags) error;
      file_to_load="";
      sread,flags(i+1),file_to_load;
      valid(i:i+1)=1;
      styc_load_buffer,file_to_load;
    }
  }
  if (anyof(valid==0)) {
    write,format=" *** ERROR: Unknow flag %s ***\n",flags(where(valid==0)(1));
    //    print_help;
  }
}

////////////////////////////////////////////////////////////
////// DISPLAY FUNCTIONS
////////////////////////////////////////////////////////////

func styc_histeq_scale(z, top=, cmin=, cmax=)
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

func styc_pltitle(title,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = pltitle_height;
  pltitle_height = long(pltitle_height*defaultdpi/83.);
  
  port= viewport();
  plth=pltitle_height;

  plt, escapechar(title), port(zcen:1:2)(1), port(4)+0.005,
    font=pltitle_font, justify="CB", height=plth;

  pltitle_height = plth_save;
}


func styc_xytitles(xtitle,ytitle,adjust,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = pltitle_height;
  pltitle_height = long(pltitle_height*defaultdpi/83.);
  
  curw=current_window();
  if (adjust==[]) {
    adjust = [0.012,0.019];
  }
  xytitles,escapechar(xtitle),escapechar(ytitle),adjust;

  pltitle_height = plth_save;
}

func escapechar(s)
{
  if (s==[]) return;
  s=streplace(s,strfind("_",s,n=20),"!_");
  s=streplace(s,strfind("^",s,n=20),"!^");
  return s;
}

func disp_1wfs_mask(nwfs)
{
  extern rtc;

  maskX = *rtc.wfs(nwfs).csX;
  maskY = *rtc.wfs(nwfs).csY;
  
  nsub = numberof(maskX);

  // FIXME ! check that refslopes are in pixels
  refsl = *rtc.wfs(nwfs).refslopes;
  if (refsl == []) return;
  
  refslX = refsl(1:nsub);
  refslY = refsl(nsub+1:);

  bX = rtc.wfs(nwfs).swX;
  bY = rtc.wfs(nwfs).swY;
  
  for (i=1;i<=nsub;i++) {
    plg,[maskY(i),maskY(i)+bY],[maskX(i),maskX(i)],marks=0,color="red";
    plg,[maskY(i),maskY(i)+bY],[maskX(i)+bX,maskX(i)+bX],marks=0,color="red";
    plg,[maskY(i),maskY(i)],[maskX(i),maskX(i)+bX],marks=0,color="red";
    plg,[maskY(i)+bY,maskY(i)+bY],[maskX(i),maskX(i)+bX],marks=0,color="red";
    //
    plmk, maskY(i)+refslY(i)+(bY-1)/2., maskX(i)+refslX(i)+(bX-1)/2., marker=2, msize=0.1, color="red";
  } 
}

func plvf(vy,vx,y,x,autoscale=,scale=,width=,hsize=,hang=,color=,type=,prop=)
/* DOCUMENT plvf,vy,vx,y,x,scale=,width=,hsize=,hang=,color=,type=,prop=
   Plots the vector field defined by (vx,vy) at positions (x,y)
   vx,vy,x,y must have the same size, but can be of arbitrary dimension.
   KEYWORDS:
   autoscale: set to 1 for the vector length to be autoscaled
   scale:     multiplicative factor applied to the autoscale results
              (for fine tweaking)
   width, color, type: same as in plg.
   hsize, hang: size and opening angle of the arrow head (default
       hsize=0.4, hang=20 degrees)
   prop:      set to zero if you want the same head size for *all* vector.
              Otherwise, the head size is proportionnal to the size of
              the vector (which results in something nicer to the eye).
   SEE ALSO: pldj
 */
{
  if (!scale) scale=1.;
  if (!width) width=2;
  if (hsize==[]) hsize=0.4;
  if (hang==[]) hang = 20;
  if (prop==[]) prop = 1;

  if (autoscale) {  
    if (prop) {
      sc=abs(vx,vy);
      if (max(sc)==0) sc=1.;
      //      else sc=sc/max(sc);
    } else {sc=1.;}

    // vector body autoscaling:
    xdif = abs(x(dif));
    w = where(xdif != 0);
    if (numberof(w)!=0) {
      minspace = min(xdif(w));
    }
    ydif = abs(y(dif));
    w = where(ydif != 0);
    if (numberof(w)!=0) {
      minspace = (minspace==[]? min(ydif(w)) : min([minspace,min(ydif(w))]) );
    }
    if (minspace==[]) minspace=1.;
    // autoscale normalization factor: max vector length / min space between location
    norm = max([vy,vx])/minspace*1.2;
    if (norm==0) norm=1.;
    vx = vx/norm*scale;
    vy = vy/norm*scale;
    //    hsize = hsize/norm*scale;
  } else {
  }
  sc = abs(vx,vy);

  pldj,(x+vx)(*),(y+vy)(*),x(*),y(*),width=width,color=color,type=type;
  x1=(x+vx)(*);  y1=(y+vy)(*);
  ang=atan(vy(*),vx(*))-(180-hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;

  ang=atan(vy,vx)-(180+hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;
}

func gen_rotation_and_noise_modes(void)
{
  extern wfs_rotation_mode,wfs_noise_mode;
  extern rtc;

  nxsub = rtc.wfs(1).sX;
  if( nxsub == 0 )
    return;
  nsubs = nxsub*nxsub;
  teldiam = 4.;
  
  xsub = (((indgen(nsubs)-1)%nxsub)+0.5)*teldiam/nxsub-teldiam/2.;
  ysub = ((indgen(nsubs)-1)/nxsub+0.5)*teldiam/nxsub-teldiam/2.;
  w = where(*rtc.wfs(1).valid);
  xsub = xsub(w);
  ysub = ysub(w);
  
  // _(x component,y component) in meters
  subpos = [xsub,ysub];
  // now in subapertures:
  subpos = subpos/teldiam*nxsub;
  // and in pixel:
  subpos = subpos*rtc.wfs(1).swX;
  
  // rotate by 1 degree:
  subposrot = transpose( mrot(1.)(+,)*subpos(,+) );

  // rotation vector field is thus defined as:
  rotv = subposrot-subpos;
  // plvf,rotv(,2),rotv(,1),subpos(,2),subpos(,1),autoscale=1;
  
  // make 1D
  rotv1d = _(rotv(,1),rotv(,2));
  // now normalize so that rotv(+)*rotv(+) = 1
  rotv1d = rotv1d/sqrt(rotv1d(+)*rotv1d(+));

  wfs_rotation_mode = rotv1d;

  // now noise mode. This should be a mode with non-zero rotational.
  // we could use the above computed rotation mode, but it has
  // unequal weight (zero in the pupil center, increasing at
  // larger pupil radius positions. Thus we can build a mode
  // with the same rotational quality but with equal radius weights:
  rotv = rotv/abs(rotv(,1),rotv(,2));
  
  // make 1D
  rotv1d = _(rotv(,1),rotv(,2));
  // now normalize so that rotv(+)*rotv(+) = 1
  rotv1d = rotv1d/sqrt(rotv1d(+)*rotv1d(+));

  wfs_noise_mode = rotv1d;
}

//---------------------------------------------------------

func mrot(ang)
/* DOCUMENT mrot(angle)
 * returns the matrix of rotation for a given angle.
 * It has to be used as follow:
 * If you want to rotate a vector of two coefficients xy=[x,y],
 * You should do rotated vector = mrot(+,)*xy(+,);
 * Angle is in degrees.
 * SEE ALSO:
 */
{
  dtor=pi/180.;
  return [[cos(ang*dtor),-sin(ang*dtor)],[sin(ang*dtor),cos(ang*dtor)]];
}

