praytop = get_env("PRAYTOP");
if (!praytop) error,"PRAYTOP is not defined!";

// load necessary files : yorick-python wrapper and styc_utils
require,praytop+"/yorick/pyk.i";
require,praytop+"/yorick/pray_utils.i";
require,praytop+"/yorick/yaodh.i";
require,praytop+"/yorick/pray_core.i";

require,praytop+"/OptimPack-1.3.2/yorick/lbfgsb.i";
require,praytop+"/OptimPack-1.3.2/yorick/OptimPack1.i";

//require,praytop+"/yorick/optimpack-mod.i";

require,praytop+"/yorick/kl.i";
require,praytop+"/yorick/widget_pray.i";

yao_fitsread = fits_read;


// start standalone version if called from shell  
arg_pray = get_argv();

if (anyof(strmatch(arg_pray,"styc"))) {
  cmd_pray = "ncpa_pray.";
 } else cmd_pray = "pray.";


if (anyof(strmatch(arg_pray,"pray.i"))) {
  cmd_pray = "";
  python_exec = praytop+"/widgets/widget_pray.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_pray  ready";
  write,"standalone version";
}

pray_wins = [26,27,28];
pray_defaultdpi = 50;    // change size of spydr graphic area
pray_itt        = 1;     // default ITT [1=lin,2=sqrt,3=square,4=log]
pray_ncolors    = 200;
pray_lut        = 0;     // default LUT index [0-41]
pray_xytitles_adjust = [0.012,0.019]; // X and Y notch axis titles in main area
pray_invertlut  = 0;
pldefault,opaque=1,marks=1;
pray_ndefoc = 0;
pray_log_itt_dex = 3;
pray_gui_realized = 0;
pray_zoom = 1;
pray_buffer_data=[];
pray_sizeimage = 100;
pray_gui = 1;
