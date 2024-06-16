#include <qi/parse_args.h>

//#ifdef HAVE_CONFIG_H
//#include "config.h"
//#endif

using namespace std;

void print_command_help (string progname) {

  cout << endl;
  cout << "Syntax is: " << progname << 
    " [-ovhn] [optional values]"  << endl;
  cout << endl;
  cout << FG_FACE1 << "   -h, --help:" << FG_NORMAL << "                     this help screen" << endl;
  cout << FG_FACE1 << "   -v, --version:" << FG_NORMAL << "                  prints version information" << endl;
  cout << FG_FACE1 << "   -n, --noise <value>:" << FG_NORMAL << "            noisify the input quadrics so that their expected size is <value>" << endl;
  cout << FG_FACE1 << "   -o, --optimize:" << FG_NORMAL << "                 Enable/Disable optimization (CPU expensive)" << endl;
  cout << FG_FACE1 << "   -aq, --affine-quadrics:" << FG_NORMAL << "         Enable/Disable print of input quadrics in affine form (in the output)" << endl;
  cout << FG_FACE1 << "   -ap, --affine-parametrizations:" << FG_NORMAL << " Enable/Disable print of parametrizations in affine form" << endl;
  cout << FG_FACE1 << "   -om, --omit:" << FG_NORMAL << "                    Enable/Disable print of components at infinity" << endl;
  cout << FG_FACE1 << "   -oe, --omit-euclidean-type" << FG_NORMAL << "      Don't show the Euclidean type of the input quadrics" << endl;
  cout << FG_FACE1 << "   -cp, --cut-params" << FG_NORMAL << "      Show the cut_params of the intersection, if any" << endl;
  cout << FG_FACE1 << "   [-l|--verbosity-level] <value>" << FG_NORMAL << "          verbosity level of the output result: 0 | 1 | 2 | 3" << endl;
  cout << "	Note: 0: prints only parametrizations and deltas  1: Brute result  2: Results with labels  3: Exhaustive results" << endl; 
  cout << FG_FACE1 << "   -q, --quiet" << FG_NORMAL << "		       invoke the qi shell without any shell messages, nor any prompt" << endl;
  cout << FG_FACE1 << "   --web" << FG_NORMAL << "			       produces the results in HTML format (QI's calculation server only)." << endl;
  cout << FG_FACE1 << "   --disable-multiline" << FG_NORMAL << "	       prevents the IO layer from producing multiline output. Useful for diff comparison during testing." << endl;
  cout << "For a quick start, just type in " << FG_FACE1 << "\"qi\"" << FG_NORMAL << ", then in the shell try the \"help\" command.\n";
  cout << endl << endl;


}

/** Check whether the libqi library has been
 *  *  compiled with the right options.
 *   *  Uses the flag "LIBQI_TESTING_COMPLIANT"
 *    *  set at configuration time */
void is_testing_compliant(void) {

	#ifdef LIBQI_TESTING_COMPLIANT
	/** Exit OK */
	exit (0);
	#else	
	/** Exit error */
	exit (1);
	#endif
}


void parseArgs (int argc, char **argv) {  // If option -o is passed, turn on optimization
  
 
  char *progname = argv[0];
 
  argc--; argv++;
  while	(argc)
    {

      if (!strcmp(*argv, "--testing-compliant"))
	      is_testing_compliant();
	    
      /* optimization flag */
      if (!strcmp(*argv, "-o") || !strcmp(*argv, "--optimize"))
	{
	  opt_optimize = true;
	  argc--; argv++;
	}

	else if (!strcmp(*argv, "-aq") || !strcmp(*argv, "--affine-quadrics"))
	{
		affineQuadrics = true;
		argc--; argv++;
	}
	
	else if (!strcmp(*argv, "-ap") || !strcmp(*argv, "--affine-parametrizations"))
	{
		affineParametrizations = true;
		argc--; argv++;
	}
	else if (!strcmp(*argv, "-oe") || !strcmp(*argv, "--omit-euclidean-type"))
	{
		showEuclideanType = false;
		argc--; argv++;
	}
	else if (!strcmp(*argv, "-cp") || !strcmp(*argv, "--cut-params"))
	{
		showCutParams = true;
		argc--; argv++;
	}
	else if (!strcmp(*argv, "-om") || !strcmp(*argv, "--omit"))
	{
		omitComponentsNotInRealSpace = true;
		argc--; argv++;
	}
	else if (!strcmp(*argv, "-l") || !strcmp(*argv, "--verbosity-level"))
	{
		argc--; argv++;
		if ( strlen(argv[0]) > 1 || (! isdigit(argv[0][0]) ) || atoi(&argv[0][0]) > 3) {
			
			cerr << FG_FACE2 << "Expecting a value in [0, 3]" << endl;
			exit (-1);
		
		}
		verbosityLevel = (VerbosityLevel)atoi(argv[0]);
		argc--; argv++;
	}
	
      /** Noisify quadrics */
      else if (!strcmp(*argv, "-n") || !strcmp(*argv, "--noise"))
	{
	  argc--; argv++;

	  // original version
	  // string_to_bigint(argv[0],noise);
	  noise = bigint(argv[0]);
	  if (noise < 2)
	    noise = 2;
	    
	  /** Initializes the random number generator */
	  // a voir pour l'initialization -> bigint.cpp 
	  //bigint::seed(time(0));
          rpl::gen.seed(static_cast<unsigned int>(time(0)));
	  argc--; argv++;
	}
      else if (!strcmp(*argv, "-q") || !strcmp(*argv, "--quiet"))
	{
		quiet = true;
		argc--; argv++;
	}
      else if (!strcmp(*argv, "--web")) {
		web = true;
		argc --; argv ++;
      }
      else if (!strcmp(*argv, "--disable-multiline")) {
		no_multiline = true;
		argc --; argv ++;
      }     
      /** Version information */
      else if (!strcmp(*argv, "-v") || !strcmp(*argv, "--version"))
	{
		//cout << LIBQI_VERSION_MAJOR << "." << LIBQI_VERSION_MINOR << "." << LIBQI_VERSION_PATCHLEVEL;
		 cout << "Version de test " <<  0.00001; 
		cout << endl;
	  exit (0);
	}
      else if (!strcmp(*argv,"-h") || !strcmp(*argv, "--help")) {
	      print_command_help (progname);
	      exit (0);
      }
      else {
	cerr << FG_FACE2 << "Unknown option: " << argv[0] << FG_NORMAL << endl;
	cerr << "Try: " << FG_FACE1 << "--help" << FG_NORMAL << " for more information." << endl << endl;
	exit (-1);
      
      }
    }

  return;

}
