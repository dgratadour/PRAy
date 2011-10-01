/*
  This collection of routines is part of the styc package
  
  Example of a yorick file associated to a widget

  Please include routines dedicated only to this widget
  Generic routines should be written in styc_utils.i file
*/

//////////////////////////////////////////////////////////
//              **********************                 //
//////    THIS IS WHERE YOU DEFINE YOUR ROUTINES   ///////
//              **********************                 //
//////////////////////////////////////////////////////////

func pray_set_starpos(type,nstars)
{
  if (type == 0) {
    pos = (random(2*nstars)-0.5)*160.;
    for (cc=1;cc<=nstars;cc++) {
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('starx%d').set_text('%s')",cc,swrite(format="%f",pos(2*cc-1)));
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('stary%d').set_text('%s')",cc,swrite(format="%f",pos(2*cc)));
    }
  }

  if (type == 1) {
    /*
    //NGS source position in mm
    mposx = [12.6,4.2,-4.2,-12.6,-12.6,-4.2,4.2,12.6,23.0,8.0,-8.0,-23.0,-23.0,-8.0,8.0,23.0,0.0,22.4,22.4,-22.4,-22.4,18.6,0.0,-18.6,0.0];
    mposy = [4.2,12.6,12.6,4.2,-4.2,-12.6,-12.6,-4.2,8.0,23.0,23.0,8.0,-8.0,-23.0,-23.0,-8.0,0.0,-22.4,22.4,-22.4,22.4,0.0,18.6,0.0,-18.6];
    // now in arcsec : 1arcsec = 621 microns
    aposx = mposx /0.621;
    aposy = mposy /0.621;
    */
    /*
    aposx = [-24.42,-37.8,-38.48,-21.2,-7.48,-13.98,-24.52,-37.86,-14.14,-7.56,-21.24,23.44,35.5,36.84,20.34,6.6,12.46,23.76,35.62,12.6,6.18,20.32,36.92,0];
    aposy = [-23.62,-36.14,-12.74,-6.56,-20.26,-37.16,24.14,37.06,38.02,21.08,7.42,24.18,36.84,13.42,7.32 ,21.06,37.94,-23.86,-36.48,-37.42,-20.92,-6.84,-13.04,0];
    */
    aposx= [-24.42,-37.8,-38.48,-21.2,-24.52,-37.86,-14.14,-7.56,-21.24,35.5,36.84,20.34,12.46,23.76,35.62,12.6,20.32,36.92,0];
    aposy = [-23.62,-36.14,-12.74,-6.56,24.14,37.06,38.02,21.08,7.42,36.84,13.42,7.32,37.94,-23.86,-36.48,-37.42,-6.84,-13.04,0];
    for (cc=1;cc<=numberof(aposx);cc++) {
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('starx%d').set_text('%s')",cc,swrite(format="%f",aposx(cc)));
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('stary%d').set_text('%s')",cc,swrite(format="%f",aposy(cc)));
    }
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_stars').set_value(%f)",19.);

  }

  if (type == 2) {
    aposx= [45.,0.,-45.,0.];
    aposy = [0.,45.,0.,-45.];
    for (cc=1;cc<=numberof(aposx);cc++) {
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('starx%d').set_text('%s')",cc,swrite(format="%f",float(aposx(cc))));
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('stary%d').set_text('%s')",cc,swrite(format="%f",float(aposy(cc))));
    }
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_stars').set_value(%f)",float(numberof(aposx)));

  }

  if (type == 3) {
    aposx= [45.,31.8198,0.,-31.8198,-45.,-31.8198,0.,31.8198];
    aposy = [0.,31.8198,45.,31.8198,0.,-31.8198,-45.,-31.8198];
    for (cc=1;cc<=numberof(aposx);cc++) {
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('starx%d').set_text('%s')",cc,swrite(format="%f",aposx(cc)));
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('stary%d').set_text('%s')",cc,swrite(format="%f",aposy(cc)));
    }
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_stars').set_value(%f)",float(numberof(aposx)));

  }
  
}


func pray_file_load(filename,images=) {
  extern pray_buffer,pray_buffer_data;
  extern pray_im,pray_name,pray_selected_display,pray_selected_star,pray_ndefoc,pray_currerr;

  if (filename != []) pray_buffer_data = fits_read(filename,fh);
  else {
    if (images != []) pray_buffer_data = images;
    else error,"nothing to do !";
  }
  
  dims = dimsof(pray_buffer_data);
  
  if (dims(2) != dims(3)) {
    pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",  "WARNING : Wrong image size. data should be square arrays.");
    return;
  }

  if (dims(1) == 3) {
    //we are dealing with classical phase diversity data
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_stars').set_value(%d)",1);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_layers').set_value(%d)",1);
    // first test if the size is a good size for ffts
    // if not add a guard band
    dims2 = fft_good(dims(2));
    if (dims2 != dims(2)) {
      mygap = long((dims2-dims(2))/2.+1);
      pray_buffer_data2 = array(0.,[3,dims2,dims2,dims(4)]);
      pray_buffer_data2(mygap:mygap+dims(2)-1,mygap:mygap+dims(2)-1,) = pray_buffer_data;
      pray_buffer_data = pray_buffer_data2;
      pray_buffer_data2 = [];
    }
    
    dims = dimsof(pray_buffer_data);  

    // checking the number of defoc positions
    pray_ndefoc = dims(4)-1;
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_defocpos').set_value(%d)",pray_ndefoc);
    if (fh != []) {
      test_styc = fits_get(fh,"STYC");
      if (test_styc != []) {
        def = fits_get(fh,"DEFOC");
        for (cc=1;cc<=pray_ndefoc/2;cc++) {
          pyk_pray,swrite(format=cmd_pray+"glade.get_widget('defoc%d').set_text('%d')",2*cc-1,cc*def);
          pyk_pray,swrite(format=cmd_pray+"glade.get_widget('defoc%d').set_text('%d')",2*cc,-cc*def);
        }
      }
    }
    pray_buffer = array(0.,[3,dims(2),dims(3),2*pray_ndefoc+2]);
    pray_buffer(,,1:pray_ndefoc+1) = pray_buffer_data;
    pray_im = array(0.,[3,dims(2),dims(2),2]);
    pray_im(,,1) = pray_buffer(,,1);
    pray_im /= max(pray_im);
    pray_im *= 10000.;
    pray_currerr = pray_buffer_data*0.;
    
    
    pray_name=["Data","Model"];
    pray_selected_display = 1;
    pray_selected_star = 0;

    for (cc=1;cc<=100;cc++)
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').remove_text(%d)",0);
    
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_sensitive(%d)",1);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').insert_text(%d,'%s')",0,"Foc");
    for (cc=1;cc<=pray_ndefoc;cc++) {
      pname = swrite(format="deFoc %d",cc);
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').insert_text(%d,'%s')",cc,pname);
    }
    
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_active(%d)",0);
    //pyk_pray,swrite(format=cmd_pray+"glade.get_widget('combo_modetype').set_active(%d)",0);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('combo_objtype').set_active(%d)",1);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('getcoeffs').set_sensitive(%d)",1);

    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
               "You have just loaded a 3-dimensional array.");
    pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","So I will run in classical mode.");

    if (pray_ndefoc >= 2)
      newstring = swrite(format="I found %d images in the cube so i assume there are %d defoc positions.",dims(4),dims(4)-1);
    else newstring = swrite(format="I found %d images in the cube so i assume there is 1 defoc position.",dims(4));
    pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
    
    if (pray_ndefoc >= 2)
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                 "Please enter the defoc values in the boxes below");
    else pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                    "Please enter the defoc value (in nm) in the box below");
    
   } else {
    if (dims(1) == 4) {
      // we are dealing with 3D phase diversity data
      dims2 = fft_good(dims(2));
      if (dims2 != dims(2)) {
        mygap = long((dims2-dims(2))/2.+1);
        pray_buffer_data2 = array(0.,[4,dims2,dims2,dims(4),dims(5)]);
        pray_buffer_data2(mygap:mygap+dims(2)-1,mygap:mygap+dims(2)-1,,) = pray_buffer_data;
        pray_buffer_data = pray_buffer_data2;
        pray_buffer_data2 = [];
      }
    
      dims = dimsof(pray_buffer_data);  

      // checking the number of defoc positions
      pray_ndefoc = dims(5)-1;
      ntarget = dims(4);
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('spin_defocpos').set_value(%d)",pray_ndefoc);
      
      pray_buffer = array(0.,[4,dims(2),dims(3),dims(4),2*pray_ndefoc+2]);
      pray_buffer(,,,1:pray_ndefoc+1) = pray_buffer_data;
      pray_im = array(0.,[3,dims(2),dims(2),2]);
      pray_im(,,1) = pray_buffer(,,1,1);
      pray_im /= max(pray_im);
      pray_im *= 10000.;
      pray_currerr = pray_buffer_data*0.;
      
      
      pray_name=["Data","Model"];
      pray_selected_display = 1;
      pray_selected_star = 1;
      
      for (cc=1;cc<=100;cc++)
        pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').remove_text(%d)",0);
      
      for (cc=1;cc<=100;cc++)
        pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_target').remove_text(%d)",0);
      
      //pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_sensitive(%d)",1);
      //pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_target').set_sensitive(%d)",1);
      
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').insert_text(%d,'%s')",0,"Foc");
      for (cc=1;cc<=pray_ndefoc;cc++) {
        pname = swrite(format="deFoc %d",cc);
        pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').insert_text(%d,'%s')",cc,pname);
      }
      
      for (cc=1;cc<=ntarget;cc++) {
        pname = swrite(format="Star # %d",cc);
        pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_target').insert_text(%d,'%s')",cc,pname);
      }
      
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_disp').set_active(%d)",0);
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_pray_target').set_active(%d)",0);
      
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('combo_modetype').set_active(%d)",0);
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('combo_objtype').set_active(%d)",1);
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('getcoeffs').set_sensitive(%d)",1);
      
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
               "You have just loaded a 4-dimensional array.");
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","So I will run in 3D mode.");
      
      if (pray_ndefoc >= 2)
        newstring = swrite(format="I found %d defoc positions.",dims(4)-1);
      else newstring = swrite(format="I found 1 defoc position.");
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",newstring);
      
      if (pray_ndefoc >= 2)
        pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                   "Please enter the defoc values in the boxes below");
      else pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                      "Please enter the defoc value (in nm) in the box below");
      
    } else {
      // wrong data format      
      pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
                 "Wrong data format. data should be either:");
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Npix x Npix x Ntarget x (Ndefoc+1) for 3D phase diversity");
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Npix x Npix x (Ndefoc+1) for classical phase diversity");
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Please load another file");
    }
  }    
  if (pray_gui) pray_setcuts;
}

func start_pray(nstars,targetx,targety,nlayers,alts,nmodes,boxsize,ndefoc,deltaFoc_nm,lambda_im,teldiam,cobs,pix_size,obj_type,obj_size,modetype,nbiter,disp,thresh,scalar,useguess,scale,diff_tt,fit_starpos,fit_object,pup_params,mir_params,script=)
{
  extern pray_buffer,pray_buffer_data;
  
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
             "Starting PRAy on the data you provided ...");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","You can abort anytime by pressing the abort button");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Please be cool and press only once ;-)");


  // parsing targetx and targety
  if (nstars > 1) {
    // dealing with 3d
    dx=strtok(targetx,"'",2*nstars);
    dy=strtok(targety,"'",2*nstars);
    xpos = ypos = array(float,nstars);
    for (rr=1;rr<=nstars;rr++) {
      sread,dx(2*rr),xpos(rr);
      sread,dy(2*rr),ypos(rr);
    }
  } else {
    xpos = [targetx];
    ypos = [targety];
  }

  // parsing nlayers and nmodes
  if (nlayers > 1) {
    dmodes=strtok(nmodes,"'",2*nlayers);
    dalts=strtok(alts,"'",2*nlayers);
    nzer = array(long,nlayers);
    alt  = array(float,nlayers);
    for (rr=1;rr<=nlayers;rr++) {
      sread,dmodes(2*rr),nzer(rr);
      sread,dalts(2*rr),alt(rr);
    }
  } else {
    nzer = [nmodes];
    alt = [alts];
  }

  // parsing deltaFoc_nm if needed
  if (ndefoc > 1) {d=strtok(deltaFoc_nm,"'",2*(ndefoc+1));
    tmp_buffer_data = array(float,ndefoc);
    for (rr=1;rr<=ndefoc;rr++) sread,d(2*rr),tmp_buffer_data(rr);
  } else {
    tmp_buffer_data = 0.;
    sread,deltaFoc_nm,tmp_buffer_data;
  }
  deltaFoc = tmp_buffer_data*2*pi/lambda_im;

  if (anyof(deltaFoc == 0)) {
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
               "---------------------------------------------------------------");
    pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Wrong defoc value entered");
    if (pray_ndefoc >= 2)
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                 "Please enter the defoc values in the boxes below");
    else pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                    "Please enter the defoc value in the box below");
    pyk_pray,cmd_pray+"y_on_loop_finished()";
    return;
  }

  if (modetype == "Zernike") tmodes = 0;
  if (modetype == "DH") tmodes = 1;
  if (modetype == "KL") tmodes = 2;
  if (modetype == "Mirror") tmodes = 3;
  if (modetype == "Canonical") tmodes = 4;
  if (modetype == "GeMS DMs") tmodes = 5;
  if (modetype == "Grid") tmodes = 6;

  // geometry init
  dims_data = dimsof(pray_buffer_data);
  if (dims_data(1) == 3) {
    images = array(float,[4,dims_data(2),dims_data(3),1,dims_data(4)]);
    images(,,1,) = pray_buffer_data;
  } else
    images = pray_buffer_data;
  
  size = dimsof(images)(2);
  // images noise variance2 init
  if (nstars > 1) var = (pray_buffer_data(1:boxsize,1:boxsize,,)(*,,)(rms,,))^2.;
  else var = (pray_buffer_data(1:boxsize,1:boxsize,)(*,)(rms,))^2.;
  //var = (pray_buffer_data(129:129+boxsize,129:129+boxsize,1)(*)(rms))^2.;

  variance2 = images*0.;
  if (nstars > 1) {
    for (kk=1;kk=nstars;kk++) {
      for (cc=1;cc<=numberof(var);cc++) {
        variance2(,,cc,kk) = clip(images(,,cc,kk),var(cc),)+var(cc);
      }
    }
  } else {
    for (cc=1;cc<=numberof(var);cc++) {
      variance2(,,cc) = clip(images(,,cc),var(cc),)+var(cc);
    }
  }
  //variance2 *= 0.5;

  if (pray_gui) pray_update_zernike_table,array(0.,nzer(1)),lambda_im,nzer(1);

  psf_core = lambda_im*1.e-9/teldiam*206265000.;//en mas
  osampl = psf_core/pix_size/2.;
  os = osampl*1.0;

  obj_nphe = sum(pray_buffer(,,1));

  if (obj_size <= 0.2) object = nphe;
  else {
    if (obj_type == "Square") {
      // square
      object=array(float,[2,size,size]);
      object(size/2,size/2)=obj_nphe/4.;
      object(size/2+1,size/2)=obj_nphe/4.;
      object(size/2,size/2+1)=obj_nphe/4.;
      object(size/2+1,size/2+1)=obj_nphe/4.;
    } else {
      object = mygauss2(size,size/2+1,size/2+1,obj_size,obj_size,max(pray_buffer(,,1)),0,0.);
      object /= sum(object);
      object *= obj_nphe;
    }
  }

  if (fit_object) fit_object = obj_size;
  
  threshold = 10.^(-thresh);
  
  for (cc=1;cc<=10;cc++)
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').remove_text(%d)",0);

  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_sensitive(%d)",1);
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",0,"Error on Foc");
  for (cc=1;cc<=pray_ndefoc;cc++) {
    pname = swrite(format="Error on deFoc %d",cc);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",cc,pname);
  }
  for (cc=numberof(deltaFoc)+1;cc<=numberof(deltaFoc)+dimsof(nzer)(2);cc++) {
    pname = swrite(format="Phase on layer %d",cc-numberof(deltaFoc));
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')", cc,pname);
  }
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_active(%d)",0);
  pray_selected_error=1;

  if (useguess) myguess = pray_param;
  
  // checking if there are a pupil map parameters
  if (pup_params != "0") {
    dx=strtok(pup_params,"'",2*7);
    pupp = array(float,7);
    for (rr=1;rr<=7;rr++) {
      sread,dx(2*rr),pupp(rr);
    }    
  } else pupp = [];

  if (mir_params != "0") {
    dx=strtok(mir_params,"'",2*4);
    mirp = array(float,4);
    for (rr=1;rr<=4;rr++) {
      sread,dx(2*rr),mirp(rr);
    }    
  } else mirp = [];

  if (scalar) res = pray(images,xpos,ypos,deltaFoc,sqrt(var),object,lambda_im,nzer=nzer,alt=alt,
                         teldiam=teldiam,cobs=cobs,osampl=os,nbiter=nbiter,tmodes=tmodes,tiptilt=tiptilt,
                         threshold=threshold,disp=disp,guess=myguess,scale=scale,diff_tt=diff_tt,
                         fit_starpos=fit_starpos,fit_object=fit_object,script=script,pup_params=pupp,mir_params=mirp);
  else {
    res = pray(images,xpos,ypos,deltaFoc,sqrt(var),object,lambda_im,nzer=nzer,alt=alt,
               teldiam=teldiam,cobs=cobs,osampl=os,nbiter=nbiter,tmodes=tmodes,tiptilt=tiptilt,
               threshold=threshold,disp=disp,variance=variance2,guess=myguess,scale=scale,diff_tt=diff_tt,
               fit_starpos=fit_starpos,fit_object=fit_object,script=script,pup_params=pupp,mir_params=mirp);
  }
}

func pray_create(nstars,targetx,targety,nlayers,alts,nmodes,ndefoc,deltaFoc_nm,lambda_im,teldiam,cobs,pix_size,obj_type,obj_size,modetype,size,snr,diff_tt,fit_starpos,fit_object,pup_params,mir_params)
{
  extern pray_data,pray_mircube;
  if (pray_gui) 
    pyk,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
               "Creating a data set with the parameters you provided");
  else
    write,"Creating a data set with the parameters you provided";
  
  // parsing deltaFoc_nm if needed
  if (ndefoc > 1) {d=strtok(deltaFoc_nm,"'",2*(ndefoc+1));
    tmp_buffer_data = array(float,ndefoc);
    for (rr=1;rr<=ndefoc;rr++) sread,d(2*rr),tmp_buffer_data(rr);
  } else {
    tmp_buffer_data = 0.;
    sread,deltaFoc_nm,tmp_buffer_data;
  }
  deltaFoc = tmp_buffer_data*2*pi/lambda_im;
  
  if (anyof(deltaFoc == 0)) {
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
               "---------------------------------------------------------------");
    pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Wrong defoc value entered");
    if (pray_ndefoc >= 2)
      pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                 "Please enter the defoc values in the boxes below and retry");
    else pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
                    "Please enter the defoc value in the box below and retry");
    return;
  }

  // parsing targetx and targety
  if (nstars > 1) {
    // dealing with 3d
    dx=strtok(targetx,"'",2*nstars);
    dy=strtok(targety,"'",2*nstars);
    xpos = ypos = array(float,nstars);
    for (rr=1;rr<=nstars;rr++) {
      sread,dx(2*rr),xpos(rr);
      sread,dy(2*rr),ypos(rr);
    }
  } else {
    xpos = [targetx];
    ypos = [targety];
  }

  // parsing nlayers and nmodes
  if (nlayers > 1) {
    dmodes=strtok(nmodes,"'",2*nlayers);
    dalts=strtok(alts,"'",2*nlayers);
    nzer = array(long,nlayers);
    alt  = array(float,nlayers);
    for (rr=1;rr<=nlayers;rr++) {
      sread,dmodes(2*rr),nzer(rr);
      sread,dalts(2*rr),alt(rr);
    }
  } else {
    nzer = [nmodes];
    alt = [alts];
  }

  if (pray_gui) pray_update_zernike_table,array(0.,100),lambda_im,nzer(1);
  
  if (modetype == "Zernike") tmodes = 0;
  if (modetype == "DH") tmodes = 1;
  if (modetype == "KL") tmodes = 2;
  if (modetype == "Mirror") tmodes = 3;
  if (modetype == "GeMS DMs") tmodes = 5;
  //if (modetype == "Grid") tmodes = 4;

  // geometry init
  psf_core = lambda_im*1.e-9/teldiam*206265000.;//en mas
  osampl = psf_core/pix_size/2.;
  os = osampl*1.;
  pupd = size/(os*2.);
  if ((long(pupd) % 2) != 0) pupd += 1;
  
  // images noise variance init
  obj_nphe = 1.e5;
  var = obj_nphe/(10.^snr);
      
  // checking if there are a pupil map parameters
  if (pup_params != "0") {
    dx=strtok(pup_params,"'",2*7);
    tmp = array(float,7);
    for (rr=1;rr<=7;rr++) {
      sread,dx(2*rr),tmp(rr);
    }    
    a = _(pupd,size/2+0.5,size/2+0.5,tmp);
    cobs = tmp(1);
    tmp = pup_model(size,a);
    pupmap = tmp/max(tmp);
  } else pupmap = [];

  pray_data             = pray_struct();
  pray_data.teldiam     = teldiam;
  pray_data.cobs        = cobs;
  pray_data.pupd        = pupd;
  pray_data.size        = size;
  pray_data.nzer        = &nzer;
  pray_data.alt         = &alt;
  pray_data.diff_tt     = diff_tt;
  pray_data.fit_starpos = fit_starpos;
  pray_data.fit_object  = fit_object;
    
  if (mir_params != "0") {
    dx=strtok(mir_params,"'",2*4);
    mirp = array(float,4);
    for (rr=1;rr<=4;rr++) {
      sread,dx(2*rr),mirp(rr);
    }    
  } else mirp = [];

  // pray init
  pray_init,xpos,ypos,tmodes=tmodes,mir_params=mirp;
  // create psfs
  if (tmodes == 2) {
    modes_coeff = _(0.,0.,gaussdev(53)/2.)(*);
  } else {
    if (tmodes == 4) {
      nzer = *prat_data.nzer;
      modes_coeff = [3.*gaussdev(nzer(1))](*);
    } else {
      modes_coeff = [7.*gaussdev(nzer(1)-1)/(indgen(nzer(1)-1)^2)](*);
    }
    //remove tt
    modes_coeff(1:2) = 0.;
    // this are the mode coefficients
    if (nlayers > 1)
      for (i=2;i<=nlayers;i++) modes_coeff=_(modes_coeff,7.*gaussdev(nzer(i))/((indgen(nzer(i))+5)^2))(*);
    modes_coeff = modes_coeff(*);
  }
  
  if (fit_starpos) coeff = _(0.1,modes_coeff(3:)); // remove tt coeffs
  else coeff = _(modes_coeff(1:2),0.1,modes_coeff(3:)); // add focus term
  // check if we want diff_tt and starpos
  if (fit_starpos) {
    starpos_coeff = gaussdev(2*nstars)*1.5; // 2 * nstars coeffs
    coeff = _(starpos_coeff,coeff);
  }
  if (diff_tt) {
    diff_coeff = random(2*ndefoc)*0.5; //2 * ndefoc coeffs
    mystr = "Adding diff TT to out of focus images";
    for (i=1;i<=2*ndefoc;i++) grow,mystr,swrite(format=" %.4f",diff_coeff(i));
    if (!pray_gui) write,mystr;
      //coeff = _(diff_coeff,coeff);
      //here we don't include it because we just run pray_coeff2psfs
  }
  if (scale) {
    scale_coeff = random(1)/5.+1.;
    //coeff = _(scale_coeff,coeff);
    //here we don't include it because we just run pray_coeff2psfs
  }
  if (fit_object) {
    obj_params = random(2)*4.;
    obj_params;
  }
  if (pray_gui) pray_update_zernike_table,coeff,lambda_im,nzer(1);
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')","Using the modes coefficients and defoc values printed below");
  
  // geometry init
  ntarget = numberof(xpos);

  if (fit_object) {
    object = mygauss2(size,size/2+1,size/2+1,1.+obj_params(1)^2,1.+obj_params(2)^2,1.,0.,0.);
    object /= sum(object);
    object *= obj_nphe;
  } else {
    if (obj_type == "Square") {
      // square
      object=array(float,[2,size,size]);
      object(size/2,size/2)=obj_nphe/4.;
      object(size/2+1,size/2)=obj_nphe/4.;
      object(size/2,size/2+1)=obj_nphe/4.;
      object(size/2+1,size/2+1)=obj_nphe/4.;
    } else {
      object = mygauss2(size,size/2+1,size/2+1,obj_size,obj_size,1.,0,0.);
      object /= sum(object);
      object *= obj_nphe;
    }
  }

  coeff_orig = coeff;
  
  if (ntarget == 1) {
    psf_tab = array(float,[3,size,size,ndefoc+1]);
    psf_tab(,,1) = pray_coeff2psfs(coeff,amp3,amp4,pupmap=pupmap);
    pray_mircube = *pray_data.mircube;
    defoc_orig = coeff(3);
    for(cc=2;cc<=ndefoc+1;cc++) {
      if (scale) coeff(3) = defoc_orig+deltaFoc(cc-1)*scale_coeff;
      else coeff(3) = defoc_orig+deltaFoc(cc-1);
      if (diff_tt) {
        coeff2 = coeff;
        coeff2(1:2) += diff_coeff(2*(cc-1)-1:2*(cc-1));
        psf_tab(,,cc) = pray_coeff2psfs(coeff2,amp3,amp4,pupmap=pupmap);
      } else psf_tab(,,cc) = pray_coeff2psfs(coeff,amp3,amp4,pupmap=pupmap);
    }
  } else {
    psf_tab = array(float,[4,size,size,ntarget,ndefoc+1]);
    psf_tab(,,,1) = pray_coeff2psfs(coeff,amp3,amp4,pupmap=pupmap);
    pray_mircube = *pray_data.mircube;
    if (fit_starpos) defoc_orig = coeff(2*ntarget+1);
    else defoc_orig = coeff(3);
    for(cc=2;cc<=ndefoc+1;cc++) {
      if (scale) {
        if (fit_starpos) coeff(2*ntarget+1) = defoc_orig+deltaFoc(cc-1)*scale_coeff;
        else coeff(3) = defoc_orig+deltaFoc(cc-1)*scale_coeff;
      } else {
        if (fit_starpos) coeff(2*ntarget+1) = defoc_orig+deltaFoc(cc-1);
        else coeff(3) = defoc_orig+deltaFoc(cc-1);
      }
      if (diff_tt) {
        coeff2 = coeff;
        coeff2(1:2) += diff_coeff(2*(cc-1)-1:2*(cc-1));
        psf_tab(,,,cc) = pray_coeff2psfs(coeff2,amp3,amp4,pupmap=pupmap);
      } else psf_tab(,,,cc) = pray_coeff2psfs(coeff,amp3,amp4,pupmap=pupmap);
    }
  }
   
  // creating images
  if (ntarget == 1) {
    images = array(float,[3,size,size,ndefoc+1]);
    for (cc=1;cc<=ndefoc+1;cc++) images(,,cc) = fft_convolve(object,eclat(psf_tab(,,cc))) + gaussdev([2,size,size])*sqrt(var);
  } else {
    images = array(float,[4,size,size,ntarget,ndefoc+1]);
    for (kk =1;kk<=ntarget;kk++) {
      for (cc=1;cc<=ndefoc+1;cc++) images(,,kk,cc) = fft_convolve(object,eclat(psf_tab(,,kk,cc))) + gaussdev([2,size,size])*sqrt(var);
    }
  }

  for (cc=1;cc<=10;cc++)
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').remove_text(%d)",0);

  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_sensitive(%d)",1);
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",0,"Error on Foc");
  for (cc=1;cc<=ndefoc;cc++) {
    pname = swrite(format="Error on deFoc %d",cc);
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')",cc,pname);
  }
  for (cc=ndefoc+1;cc<=ndefoc+dimsof(pray_mircube)(4);cc++) {
    pname = swrite(format="Phase on layer %d",cc-numberof(deltaFoc));
    pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').insert_text(%d,'%s')", cc,pname);
  }
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('winselect_error_disp').set_active(%d)",0);
  pray_selected_error=1;

  pray_file_load,images = images;

  return coeff_orig;
}


//////////////////////////////////////////////////////////
// A basic set of display functions
//////////////////////////////////////////////////////////
func pray_disp_select(ndisp,nstar)
{
  extern pray_selected_display,pray_selected_star,pray_im,pray_name,pray_ndefoc,pray_buffer;

  if (ndisp != []) pray_selected_display = ndisp + 1;
  if (nstar != []) pray_selected_star = nstar + 1;

  if (pray_selected_star == 0) {
    pray_im(,,1) = pray_buffer(,,1+(pray_selected_display-1));
    pray_im(,,2) = pray_buffer(,,1+pray_ndefoc+1+(pray_selected_display-1));
  } else {
    pray_im(,,1) = pray_buffer(,,pray_selected_star,1+(pray_selected_display-1));
    pray_im(,,2) = pray_buffer(,,pray_selected_star,1+pray_ndefoc+1+(pray_selected_display-1));
  }
  pray_name=["Data","Model"];
  
  pray_im /= clip(max(pray_im),1.e-6,);
  pray_im *= 10000.;
  
  pray_setcuts;
  //pray_disp_error,pray_selected_display-1;
}

func pray_disp_error(ndisp)
{
  extern pray_selected_error,pray_currerr,pray_nlayers,pray_mircube,pray_ndefoc;
  extern pray_zoom,dispok;

  if (ndisp != []) pray_selected_error = ndisp + 1;
  if (pray_selected_error == []) return;
  
  window,pray_wins(3);
  fma;
  
  if (pray_selected_error > pray_ndefoc + 1) {
    //pray_mircube = (*pray_data.mircube)*(*(pray_data.ipupil));
    // then we want to display the phase per layer
    nlayer = pray_selected_error - (pray_ndefoc + 1);
    tmpdiam = (pray_data.pupd + nlayer-1)/2;
    tmpdiam2 = pray_data.size/2;
    pli,pray_mircube(tmpdiam2-tmpdiam:tmpdiam2+tmpdiam+1,tmpdiam2-tmpdiam:tmpdiam2+tmpdiam+1,nlayer);
    tmptlt = swrite(format="Phase on layer %d",nlayer);
  } else {
    // then we want to display the error for the selected star
    if (pray_zoom != 1) {
      npix = dimsof(pray_currerr)(2);
      newpix = long(npix/pray_zoom);
      mygap = long((npix-newpix)/2.+1);
      if (pray_selected_star != 0) pli,pray_currerr(mygap:mygap+newpix-1,mygap:mygap+newpix-1,pray_selected_star,pray_selected_error);
      else pli,pray_currerr(mygap:mygap+newpix-1,mygap:mygap+newpix-1,pray_selected_error); 
    } else {
      if (pray_selected_star != 0) pli,pray_currerr(,,pray_selected_star,pray_selected_error);
      else pli,pray_currerr(,,pray_selected_error);
    }
    
    if (pray_selected_error == 1) tmptlt = "Error on Foc";
    else tmptlt = swrite(format="Error on deFoc %d",pray_selected_error-1);
  }
  pray_pltitle,tmptlt,defaultdpi=pray_defaultdpi;
  colorbar,adjust=-0.024,levs=10;
}

func pray_progressbar_frac(frac)
{
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('progressbar_pray').set_fraction(%f)",float(frac));  
}

func pray_progressbar_text(text)
{
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('progressbar_pray').set_text('%s')",text);  
}

func pray_dispgui_update
{
  //pyk_pray,swrite(format="y_pray_set_lut(%d)",pray_lut);
  pyk_pray,swrite(format=cmd_pray+"pray_set_itt(%d)",clip(long(pray_itt-1),0,4));

}

func pray_win_init(pid1,pid2,pid3,redisp=)
{
  extern pray_xid1,pray_xid2,pray_xid3;
  extern pray_defoc,pray_defaultdpi,pray_lut,pray_wins,pray_gui_realized;
  
  pray_xid1=pid1; pray_xid2=pid2;pray_xid3=pid3;

  if (catch(0x08)) {
    pray_gui_realized = 1;
  }
  
  if (!pray_gui_realized) {
    stylename = praytop+"/gs/pray.gs";
 
    window,pray_wins(1),dpi=pray_defaultdpi,width=0,height=0,xpos=-2,ypos=-2,style=stylename,parent=pid1;
    limits,square=1;
    palette,"gray.gp"; // need this if loadct is used!?
    
    window,pray_wins(2),dpi=pray_defaultdpi,width=0,height=0,xpos=-2,ypos=-2,style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp"; // need this if loadct is used!?
    
    window,pray_wins(3),dpi=pray_defaultdpi,width=0,height=0,xpos=-2,ypos=-2,style=stylename,parent=pid3;
    limits,square=1;
    palette,"gray.gp"; // need this if loadct is used!?
    
    pray_set_lut,pray_lut;
    
    if (redisp) return;
    
    if (pray_ndefoc>0) {
      pray_setcuts;
    }
  }
  pray_gui_realized = 1;
  pyk_pray,swrite(format=cmd_pray+"glade.get_widget('textview2').get_buffer().set_text('%s')",
             "Welcome to PRAy v0.9 : the Phase Retrival Algorithm in Yorick");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "---------------------------------------------------------------");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "Please load some data using the dedicated button");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "You can load either a N x N x Ndefoc array");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "                    or a N x N x Ntargets x Ndefoc array");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "If each case you have to fill the defoc box(es) below");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "In the 3d case you also have to provide targets positions in the FoV");
  pyk_pray,swrite(format=cmd_pray+"y_add_comment_txt_nl('%s')",
             "Enjoy ...");
  
  
}

func pray_set_lut(_lut)
// change the LookUp Table and Intensity Transfer Table
{
  require,"idl-colors.i";
  extern rlut,glut,blut;
  extern pray_invertlut,pray_lut,pray_itt,pray_ncolors,pray_log_itt_dex;
  local r,g,b;

  //  if ((!lut)||(lut==itp_lut)) return; // nothing to do.

  window,pray_wins(1);
  if (_lut!=[]) pray_lut = _lut;
  
  if (_lut!=[]) {  // then read and set new lut
    if (_lut==0) palette,"earth.gp";
    else loadct,_lut;
    palette,query=1,rlut,glut,blut;  // store
  }

  // invert?
  if (pray_invertlut) {
    r=rlut(::-1); g=glut(::-1); b=blut(::-1);
  } else {
    r=rlut; g=glut; b=blut;
  }

  // itt:
  if (pray_itt<=1) { // linear
    ind = span(0.,1.,pray_ncolors);
  } else if (pray_itt==2) { // sqrt
    ind = sqrt(span(0.,1.,pray_ncolors));
  } else if (pray_itt==3) { // square
    ind = (span(0.,1.,pray_ncolors))^2.;
  } else if (pray_itt==4) { // log
    ind = log10(span(10.^(-pray_log_itt_dex),1.,pray_ncolors)); // 8 dex
    ind -= min(ind);
    ind /= max(ind);
  } else if (pray_itt>=5) { // histeq
    ind = span(0.,1.,pray_ncolors);
  }
  ind = long(round(ind*(pray_ncolors-1)+1));
  r = r(ind); g = g(ind); b = b(ind);

  // and finally, load the palette:
  for (i=1;i<=3;i++) {
    window,pray_wins(i);
    palette,r,g,b;
  }
  pray_disp;
}

func pray_disp(void)
// pli display, main window
{
  extern pray_ndefoc,pray_wins;
  extern pray_im,pray_imd,pray_imdnum,pray_cmin,pray_cmax,pray_name;
  extern pray_zoom;

  if (pray_ndefoc == 0) return;
  //  write,format="%s ","*";

  pray_imd = pray_im*0.;
  
  for (i=1;i<=2;i++) {
    if (pray_itt==5) {
      if (nallof(pray_imdnum==[2,pray_itt]))     \
        pray_imd(,,i) = pray_histeq_scale(pray_im(,,i));
    } else pray_imd(,,i) = bytscl(pray_im(,,i),cmin=pray_cmin,cmax=pray_cmax);

    window,pray_wins(i);
    fma;
    if (pray_zoom != 1) {
      npix = dimsof(pray_im)(2);
      newpix = long(npix/pray_zoom);
      mygap = long((npix-newpix)/2.+1);
      pli,pray_imd(mygap:mygap+newpix-1,mygap:mygap+newpix-1,i);
    } else pli,pray_imd(,,i);

    axtit = "pixels";
    
    pray_pltitle,pray_name(i),defaultdpi=pray_defaultdpi;
    pray_xytitles,axtit,axtit,defaultdpi=pray_defaultdpi;
    colorbar,adjust=-0.024,levs=10;
  }
  
  pray_imdnum = [2,pray_itt];
  pray_dispgui_update;
}


func pray_disp_cpc(e,all=)
{
  extern pray_cmin,pray_cmax,pray_ndefoc;

  if (pray_ndefoc==0) return;

  if (max(pray_im(,,2)) > 0.) tmp = 2;
  else tmp = 1;
                                
  if (all) eq_nocopy,subim,pray_im(,,tmp);
  else subim=pray_get_subim(x1,x2,y1,y2,tmp);
  
  if (subim==[]) subim = pray_im(,,tmp); // init, not defined
  
  if (e==0) {
    cmin = min(subim);
    cmax = max(subim);
  } else {
    cmin = min(pray_cpc(subim,0.1,0.999));
    cmax = max(pray_cpc(subim,0.1,0.999));
  }
  
  pray_cmin = cmin;
  pray_cmax = cmax;
  pyk_pray,swrite(format=cmd_pray+"pray_set_cmincmax(%f,%f,%f)",float(cmin),float(cmax),float(cmax-cmin)/100.);
}

func pray_cpc(im,fmin,fmax)
/* DOCUMENT func cpc(im,fmin,fmax)
   return clipped image at
   from cmin=fraction fmin of pixel intensity
   to   cmax=fraction fmax of pixel intensity
   0 <= fmin < fmax <= 1
   example: pli,cpc(im)
   SEE ALSO:
 */
{
  s = sedgesort(im(*));
  n = numberof(im);
  if (fmin==[]) fmin = 0.10;
  if (fmax==[]) fmax = 0.995;
  if ((fmin<0)||(fmax<0)) error,"fmin and fmax should be > 0";
  if ((fmin>1)||(fmax>1)) error,"fmin and fmax should be < 1";
  if (fmin>fmax) error,"fmin should be < fmax";
  x1=s(long(round(n*fmin)));
  x2=s(long(round(n*fmax)));
  return clip(im,x1,x2);
}

func pray_set_cmin(pycmin)
{
  extern pray_cmin,pray_cmax;
  //  write,format="pycmin=%f, cmin=%f\n",pycmin,cmin;
  // because the precision was cut by the yorick -> python -> yoric
  // transfer, we have to allow for some slack
  if ((pray_cmax-pray_cmin)!=0) if ((abs(pycmin-pray_cmin)/(pray_cmax-pray_cmin))<1e-3) return;
  if (pycmin>pray_cmax) {
    write,"cmin > cmax, ignoring";
    return;
  }
  pray_cmin = pycmin;
  pray_disp;
}

func pray_set_cmax(pycmax)
{
  extern pray_cmax,pray_cmin;
  //  write,format="pycmax=%f, cmax=%f\n",pycmax,cmax;
  // because the precision was cut by the yorick -> python -> yoric
  // transfer, we have to allow for some slack
  if ((pray_cmax-pray_cmin)!=0) if ((abs(pycmax-pray_cmax)/(pray_cmax-pray_cmin))<1e-3) return;
  if (pycmax<pray_cmin) {
    write,"cmax < cmin, ignoring";
    return;
  }
  pray_cmax = pycmax;
  pray_disp;
}

func pray_get_subim(&x1,&x2,&y1,&y2,nwin)
{
  extern pray_ndefoc;
  
  if (pray_ndefoc==0) return;
  curw = current_window();
  window,pray_wins(nwin);
  lim=limits();
  dims = dimsof(pray_im(,,nwin));
  x1=long(floor(clip(lim(1),1,dims(2))));
  x2=long(ceil(clip(lim(2),1,dims(2))));
  y1=long(floor(clip(lim(3),1,dims(3))));
  y2=long(ceil(clip(lim(4),1,dims(3))));
  if ( (x1==x2) || (y1==y2) ) {
    //    itp_pyk_status_push,"WARNING: (get_subim) Nothing to show";
    //    write,"WARNING: (get_subim) Nothing to show";
    return;
  }
  window,curw;
  return pray_im(x1:x2,y1:y2,nwin)
}

func pray_unzoom
{
  extern pray_zoom;
  
  for (i=1;i<=2;i++) {
    window,pray_wins(i);
    unzoom;
    limits;
  }
  pray_zoom = 1;
}

func pray_setcuts
{
  pray_disp_cpc,all=1;
  pray_disp;
}

func pray_send_to_spydr
{
  spydr,pray_buffer;
}

func pyk_pray(cmd)
{

  if (pray_gui != 0) pyk,cmd;
}
//////////////////////////////////////////////////////////
//              **********************                 //
////////         END OF ROUTINES DEFINITION      ////////
//              **********************                 //
//////////////////////////////////////////////////////////

func pray_quit
{
  extern always_busy;
  
  if (anyof(strmatch(arg_stats,"widget_slodar.i"))) {
    quit;
  } else {
    always_busy = 0;
    worker_ready;
    quit;
  }
}


//22 24;
//23 27
