praytop = get_env("PRAYTOP");
if (!praytop) error,"PRAYTOP is not defined!";

require,praytop+"/yorick/pray.i";

write,"pray env loaded !";
pray_gui = 0;

/*
func pray_create(nstars,targetx,targety,nlayers,alts,nmodes,ndefoc,deltaFoc_nm,lambda_im,teldiam,cobs,pix_size,obj_type,obj_size,modetype,size,snr,diff_tt,fit_starpos,fit_object)  
 */


func test2d(snr,ndefoc,defoc)
{
  extern mcoeff,ref_phase;

  mcoeff = pray_create(1,0.,0.,1,0.,55,ndefoc,defoc,1550.,7.9,0.,15.,"Gaussian",1.,"Zernike",64,snr,0,0,0);
  ref_phase = pray_mircube;
  write,"starting pray loop";
  start_pray,1,0.,0.,1,0.,55,30,ndefoc,defoc,1550.,7.9,0.,15.,"Gaussian",1.,"Zernike",500,1,3,1,0,0,0,0,0,script=1;
  //(pray_mircube(,,1)-ref_phase(,,1))(*)(rms)*1600./pi/2.
  valids = where(*pray_data.ipupil);
  write,format="Initial phase rms (nm)  : %f\n",ref_phase(,,1)(valids)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) : %f\n",(ref_phase(,,1)-pray_mircube(,,1))(valids)(*)(rms)*1.e3;
  return (ref_phase(,,1)-pray_mircube(,,1))(valids)(*)(rms)/ref_phase(,,1)(valids)(*)(rms);
}

func test2dscale(snr)
{
  extern mcoeff,ref_phase;

  mcoeff = pray_create(1,0.,0.,1,0.,55,1,"400",1550.,7.9,0.,15.,"Gaussian",1.,"Zernike",64,snr,0,0,1);
  ref_phase = pray_mircube;
  write,"starting pray loop";
  start_pray,1,0.,0.,1,0.,55,30,1,"400",1550.,7.9,0.,15.,"Gaussian",1.,"Zernike",500,1,3,1,0,0,0,0,1,script=1;
  //(pray_mircube(,,1)-ref_phase(,,1))(*)(rms)*1600./pi/2.
  valids = where(*pray_data.ipupil);
  write,format="Initial phase rms (nm)  : %f\n",ref_phase(,,1)(valids)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) : %f\n",(ref_phase(,,1)-pray_mircube(,,1))(valids)(*)(rms)*1.e3;
  return (ref_phase(,,1)-pray_mircube(,,1))(valids)(*)(rms)/ref_phase(,,1)(valids)(*)(rms);
}

func script2d(void)
{

  tabres = array(0.,10,100,5);
  tabdefoc = ["100","200","300","400","500","600","700","800","900","1000"];
  for (i=1;i<=10;i++) {
    for (j=1;j<=5;j++) {
      for (k=1;k<=100;k++) tabres(i,k,j) = test2d(j+1,1,tabdefoc(j));
    }
  }

}



func test3d(snr)
{
  extern mcoeff,ref_phase;

  mcoeff=pray_create(19,"['-24.420000', '-37.800000', '-38.480000', '-21.200000', '-24.520000', '-37.860000', '-14.140000', '-7.560000', '-21.240000', '35.500000', '36.840000', '20.340000', '12.460000', '23.760000', '35.620000', '12.600000', '20.320000', '36.920000', '0.000000']", "['-23.620000', '-36.140000', '-12.740000', '-6.560000', '24.140000', '37.060000', '38.020000', '21.080000', '7.420000', '36.840000', '13.420000', '7.320000', '37.940000', '-23.860000', '-36.480000', '-37.420000', '-6.840000', '-13.040000', '0.000000']",3,"['0', '4500', '9000']","['296', '416', '208']",2,"['-300', '300']",1550.,7.9,0.,15.,"Gaussian",1.,"Zernike",64,snr,0,0,0);
  
  ref_phase = pray_mircube;
  
  start_pray,19,"['-24.420000', '-37.800000', '-38.480000', '-21.200000', '-24.520000', '-37.860000', '-14.140000', '-7.560000', '-21.240000', '35.500000', '36.840000', '20.340000', '12.460000', '23.760000', '35.620000', '12.600000', '20.320000', '36.920000', '0.000000']", "['-23.620000', '-36.140000', '-12.740000', '-6.560000', '24.140000', '37.060000', '38.020000', '21.080000', '7.420000', '36.840000', '13.420000', '7.320000', '37.940000', '-23.860000', '-36.480000', '-37.420000', '-6.840000', '-13.040000', '0.000000']", 3, "['0', '4500', '9000']", "['296', '416', '208']", 30, 2, "['-300', '300']", 1550.,7.9, 0., 15., "Gaussian",1., "Zernike", 500, 1, 3, 1, 0, 0, 0, 0, 0,script=1;

  for (i=1;i<=3;i++) {
    write,format="Initial phase rms (nm)  on layer %d : %f\n",i,ref_phase(,,1)(*)(rms)*1.e3;
    write,format="Residual phase rms (nm) on layer %d : %f\n",i,(ref_phase(,,1)-pray_mircube(,,1))(*)(rms)*1.e3;
  }
}

/*
to analyze the perf we need to have reference cases.
relatively low altitude layers so that we have a good overlaping

4 stars +- 45" appart
this means h > 1/2 diam / theta
(so that the footprint of the metapup is half the diam larger)
in the case of 45" and 8m tel this leads to :
4. /(45 * 4.848e-6) = 18335.2 m
if we take half of this we have a good overlapping

so we take 4 stars :
"['-45.','0.','45.','0.']"
"['0.','-45.','0.','45.']"
and 2 layers : 
"['0','9000']"
one defoc : "400"




*/



func test3dsimple(snr)
{
  extern mcoeff,ref_phase;

  mcoeff=pray_create(4,"['-45.','0.','45.','0.']","['0.','-45.','0.','45.']",
                     2,"['0', '4500']","['100', '100']",
                     1,"400",
                     1550.,7.9,0.,15.,"Gaussian",.5,"Zernike",64,snr,0,0,0);
  
  ref_phase = pray_mircube;
  
  start_pray,4,"['-45.','0.','45.','0.']","['0.','-45.','0.','45.']",
  2, "['0', '4500']", "['100', '100']", 15,
  1, "400",
  1550.,7.9, 0., 15., "Gaussian",.5, "Zernike", 500, 1, 3, 1, 0, 0, 0, 0, 0,script=1;

  tmp = *pray_data.def;
  write,format="Initial phase rms (nm)  on layer %d : %f\n",1,ref_phase(,,1)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) on layer %d : %f\n",1,(ref_phase(,,1)-tmp(,,1:100)(,,+)*pray_param(1:100)(+))(*)(rms)*1.e3;

  write,format="Initial phase rms (nm)  on layer %d : %f\n",2,ref_phase(,,2)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) on layer %d : %f\n",2,(ref_phase(,,2)-tmp(,,101:)(,,+)*pray_param(101:)(+))(*)(rms)*1.e3;
  
}

// 8 stars ...

func test3dsimple2(snr,defoc)
{
  extern mcoeff,ref_phase;

  test_rms = 400.;
  while (test_rms > 300) {
    mcoeff=pray_create(8,"['45.000000', '31.819800', '0.000000', '-31.819800', '-45.000000', '-31.819800', '0.000000', '31.819800']", "['0.000000', '31.819800', '45.000000', '31.819800', '0.000000', '-31.819800', '-45.000000', '-31.819800']", 2, "['0', '9000']", "['100', '100']", 1, defoc, 1550.0, 7.9, 0.12, 15.0, "Gaussian", 1.0, "Zernike", 64, snr, 0, 0, 0, "0", "0");

  
    ref_phase = pray_mircube;
    test_rms = ref_phase(,,1)(*)(rms)*1.e3;
  }


  
  start_pray, 8, "['45.000000', '31.819800', '0.000000', '-31.819800', '-45.000000', '-31.819800', '0.000000', '31.819800']", "['0.000000', '31.819800', '45.000000', '31.819800', '0.000000', '-31.819800', '-45.000000', '-31.819800']", 2, "['0', '9000']", "['100', '100']", 30, 1, defoc, 1550.000000, 7.900000, 0.120000, 15.000000, "Gaussian", 1.000000, "Zernike", 500, 1, 3, 1, 0, 0, 0, 0, 0, "0", "0",script=1;

  tmp = *pray_data.def;
  
  write,format="Initial phase rms (nm)  on layer %d : %f\n",1,ref_phase(,,1)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) on layer %d : %f\n",1,(ref_phase(,,1)-tmp(,,1:100)(,,+)*pray_param(1:100)(+))(*)(rms)*1.e3;

  write,format="Initial phase rms (nm)  on layer %d : %f\n",2,ref_phase(,,2)(*)(rms)*1.e3;
  write,format="Residual phase rms (nm) on layer %d : %f\n",2,(ref_phase(,,2)-tmp(,,101:)(,,+)*pray_param(101:)(+))(*)(rms)*1.e3;

  return [(ref_phase(,,1)-tmp(,,1:100)(,,+)*pray_param(1:100)(+))(*)(rms)/ref_phase(,,1)(*)(rms),
          (ref_phase(,,2)-tmp(,,101:)(,,+)*pray_param(101:)(+))(*)(rms)/ref_phase(,,2)(*)(rms)]
  
}

func test3dsimple_full(void)
{

  ndefoc = 4;
  nsnr = 5;
  tabdefoc = ["200","400","600","800"];
  tabres = array(0.,nsnr,ndefoc,100,2);
  for (i=1;i<=nsnr;i++) {
    for (j=1;j<=ndefoc;j++) {
      for (k=1;k<=100;k++) tabres(i,j,k,) = test3dsimple2(i+1,tabdefoc(j));
    }
  }
  return tabres;
}
// quadratic modes filtering
//
