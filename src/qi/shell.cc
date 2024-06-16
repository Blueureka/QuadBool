#include <qi/shell.h>
#include <sstream>


using namespace std;
using namespace QI;
using namespace rpl;

static string prompt = FG_BOLD  "qi# "  FG_NORMAL ;
static ShellVar shell_vars[MAX_VARIABLES];
static short current_variable = 0;
static bool useLaTeX = false;

void print_shell_help () {

  cout << endl;
  cout << "Available commands are:" << endl << endl;
  cout << FG_FACE1 << "  quit" << FG_NORMAL << endl;
  cout << FG_FACE1 << "  intersect " << FG_NORMAL << "<var1> <var2>    : intersects two quadrics <var1> and <var2>" << endl;
  cout << FG_FACE1 << "  noise " << FG_NORMAL << "<value>              : defines a noise value to perturb the intersection calculation" << endl;
  cout << FG_FACE1 << "  affine-quad" << FG_NORMAL << "                : Enable/Disable print of input quadrics in affine form (in the output)" << endl;
  cout << FG_FACE1 << "  affine-param" << FG_NORMAL << "               : Enable/Disable print of parametrizations in affine form" << endl;
  cout << FG_FACE1 << "  euclidean-type" << FG_NORMAL << "             : Enable/Disable display of the Euclidean type of the input quadrics" << endl;
  cout << FG_FACE1 << "  omit" << FG_NORMAL << "                       : Enable/Disable print of components at infinity" << endl;
  cout << FG_FACE1 << "  optimize" << FG_NORMAL << "                   : Enable/Disable optimization (CPU expensive)" << endl;
  cout << FG_FACE1 << "  level" << FG_NORMAL << "                      : Defines the verbosity level [0, 2]" << endl;
  cout << "	Note: 0: prints only parametrizations and deltas 1: Brute result  2: Results with labels  3: Exhaustive results" << endl;
  cout << FG_FACE1 << "  !" << FG_NORMAL << "                          : Repeat last command" << endl;
  cout << endl;
  cout << "To set a variable, type, for instance:" << endl;
  cout << endl <<"  my_quadric=x^2-x+3*y*z-50" << endl;
  cout << endl;
  cout << " Note: To define a quadric, you can use the following representations:" << endl;
  cout << "          - Polynomial                ex: x^2-x+3*y*z-50" << endl;
  cout << "          - Homogeneous polynomial    ex: x^2-x*w+3*y*z-50*w^2" << endl;
  cout << "          - Matricial                 ex: [1 0 0 -1 0 3 0 0 0 -50]" << endl << endl;
  cout << " Matricial coefficients are stored as follow:" << endl;
  cout << " [x^2 xy xz xw y^2 yz yw z^2 zw w^2]" << endl << endl;
  
}

short index_of(char *varname) {

  for (short i = 0; i < MAX_VARIABLES; i++) 
    { 
      if ( ! strcmp(shell_vars[i].name, varname) )
	return i;
    }

  return -1;

}
 
void set_value_of(char *varname, char *value) {

  short indice = -1;
  QIParser parser;

  parser.setQuadricDesc (string(value));

  try {

    parser.parse();
    
    if ( (indice = index_of(varname)) < 0 ) {
      if ( ! quiet)
      	cout << "Adding a new variable named: " << FG_FACE3 << varname << FG_NORMAL << endl;
      indice = current_variable;
      current_variable++;
      current_variable %= MAX_VARIABLES;
      strcpy (shell_vars[indice].name, varname);
    }
    
    if ( ! quiet) {
      cout << "Affecting variable: " << varname << " with value: " << FG_FACE3 << value << FG_NORMAL << endl;
      cout << "Vectorial desc:     " << parser.getVectorialDesc() << endl;
    }
    
    shell_vars[indice].value = parser.getMatricialDesc();
    


  }catch (string parse_error) { cout << parse_error << endl; }

}

void noisify (bigint_matrix *q1, bigint_matrix *q2) {

  if ( ! quiet)
   cout << FG_FACE1 << "Noisifying ...";

  bigint rand_size = cubic_root(noise);

  bigint_matrix p = rand_mat(rand_size);
  while (det(p) == 0)
    p = rand_mat(rand_size);

  bigint r1,r2,r3,r4;

  r1 = rand_bigint(rand_size);
  r2 = rand_bigint(rand_size);
  r3 = rand_bigint(rand_size);
  r4 = rand_bigint(rand_size);

  while (r1*r4-r2*r3 == 0)
    {
      r1 = rand_bigint(rand_size);
      r2 = rand_bigint(rand_size);
      r3 = rand_bigint(rand_size);
      r4 = rand_bigint(rand_size);
    }

  bigint_matrix q1_r = prod(base_matrix<bigint> (prod(trans(p), (*q1))), p);
  bigint_matrix q2_r = prod(base_matrix<bigint> (prod(trans(p), (*q2))), p);
		  
  *q1 = r1*q1_r+r2*q2_r;
  *q2 = r3*q1_r+r4*q2_r;

  if ( ! quiet)
  	cout << "Done" << FG_NORMAL << endl;

}

void shell_main (void) {

   #ifndef NDEBUG
     std::cout << "entering  shell_main" << std::endl;
   #endif 
   
   if ( ! quiet ) {
    //cout << endl << FG_FACE3 << "QI shell" << FG_NORMAL << endl << "(based on libqi-" << LIBQI_VERSION_MAJOR << "." << LIBQI_VERSION_MINOR << "." << xstringify(LIBQI_VERSION_PATCHLEVEL);
    cout << endl << FG_FACE3 << "QI shell" << FG_NORMAL << endl << "(based on libqi-" << 0.00001; 
    cout << ")" << endl << endl;
    cout << "Type " << FG_FACE1 << "\"help\"" << FG_NORMAL << " to know what to do" << endl;
   }
  
   char command[MAX_COMMAND_SIZE];
   char last[MAX_COMMAND_SIZE];
   
   for ( ; ; ) {

   if ( ! quiet )
     cout << prompt;
    
    memset(command, 0, MAX_COMMAND_SIZE*sizeof(char));
//    char * ret = fgets(command, MAX_COMMAND_SIZE, stdin);
fgets(command, MAX_COMMAND_SIZE, stdin);
    // on ne regarde pas la valeur de retour car 
    // le parser doit joue la dessus pour s'arreter
    /** Remove trailing \n */
    command[strlen(command)-1] = '\0';

    /** Ignore comment lines, such as the magic shebang line of Unix shells. */
    if ( command[0] == '#' ) continue;
    
    if ( ! strcmp(command, "!") ) strncpy (command, last, MAX_COMMAND_SIZE);
    else strncpy(last, command, MAX_COMMAND_SIZE);
    
    if       ( ! strcmp(command, "help") ) print_shell_help();
    else if  ( ! strcmp(command, "quit") || feof(stdin) ) exit(0);
    else if  (   strstr(command, "intersect") ) {
      
      char *token = command;
      char var1[16];
      char var2[16];
      
      var1[0]='\0';
      var2[0]='\0';

      strtok (command, " ");
      while ( (token = strtok(NULL, " ")) ) {

	if      ( var1[0] == '\0' ) strcpy(var1, token);
	else if ( var2[0] == '\0' ) strcpy(var2, token);
	else {
	  if ( ! quiet)
	   cout << "Ignoring extra arguments: " << token << endl;
	}
      }

      if ( var1[0] == '\0' || var2[0] == '\0') {
	cout << FG_FACE2 << "intersect: Arguments missing." << FG_NORMAL << endl;
	continue;
      }

      bigint_matrix q1;
      bigint_matrix q2;
      string noise_str = "";
      short i;

      if ( (i = index_of (var1)) < 0 ) {
	cerr << FG_FACE2 << "Undefined variable: " << var1 << FG_NORMAL << endl;
	continue;
      }
      q1 = shell_vars[i].value;

      if ( (i = index_of (var2)) < 0 ) {
	cerr << FG_FACE2 << "Undefined variable: " << var2 << FG_NORMAL << endl;
	continue;
      }
      q2 = shell_vars[i].value;

      if (noise != 0) {
	char *tmp;
	tmp = new char[MAX_COMMAND_SIZE];
	noisify (&q1, &q2);
        // original version
        tmp = bigint_to_string (noise);
//	noise = bigint(tmp);
	noise_str = string(tmp);
        delete [] tmp;
      }
      
      if ( ! quiet ) {
       cout << FG_FACE1 << "Intersecting " << FG_FACE3 << var1 << " (";
       print_quadric(q1, cout);
       cout << ")" << FG_NORMAL << " and " << FG_FACE3 << var2 << " (";
       print_quadric(q2, cout);
       cout << ")" << FG_NORMAL << endl << endl;
      }
      
      quad_inter <bigint> qi_ic = intersection(q1,q2,opt_optimize,cout);
	
      QIOutputter outputter;
      QIWriter *writer = NULL;

      if ( web )
      	writer = new QIHTMLWriter ();

      if ( writer == NULL )
	writer = new QIConsoleWriter ();
      
      if (affineQuadrics)
      	outputter.setAffineInputQuadrics();
      if (affineParametrizations)
        outputter.setAffineParametrizations();
      if (omitComponentsNotInRealSpace)
        outputter.setOmitImaginaryParametrizations();
      if (showEuclideanType)
	outputter.showInputEuclideanType();
      if (showCutParams)
	outputter.showCutParams();
      
      outputter.useLaTeX(useLaTeX);
      outputter.setVerbosityLevel(verbosityLevel);

      if (no_multiline)
	      outputter.disableMultiline ();
      
      /** Transform kernel results into strings */
      outputter.output (qi_ic, q1, q2);
 
      /** Console rendering */
      writer->setOutputInformation(outputter.getOutput());
      writer->setVerbosityLevel(verbosityLevel);
      writer->setNoiseValue (noise_str);
      writer->write ();/** Produces the result to /dev/stdout by default */
      
      delete writer;
      
    }
    else if ( strstr (command, "=") ) {

      char *varname = NULL;
      char *value   = NULL;
      
      varname = strtok (command, "=");
      value   = strtok (NULL, "=");
      
      if ( ! value || strlen(value) == 0)
	cout << FG_FACE2 << "You must provide a value" << FG_NORMAL <<endl;
      else
	set_value_of (varname, value);

      
    }
    else if ( strstr (command, "noise") ) {

      char *value = NULL;
	
      if ( ! quiet)
       cout << "Old noise value: " << FG_FACE3 << noise << FG_NORMAL << endl;
      strtok (command, " ");

      value = strtok(NULL, " ");

      if ( ! value || strlen(value) == 0) {
	cout << FG_FACE2 << "noise: You must provide a value" << FG_NORMAL << endl;
	continue;
      }
      else
	// original version
	//string_to_bigint (value, noise);
        noise = bigint(value); 

      if ( ! quiet)
       cout << "New noise value: " << FG_FACE3 << noise << FG_NORMAL << endl;

    }
    else if ( strstr (command, "level") ) {

	    char *value = NULL;
	    
	    if ( ! quiet)
	     cout << "Old verbosity level: " << FG_FACE3 << verbosityLevel << FG_NORMAL << endl;
	    
	    strtok (command, " ");

	    value = strtok(NULL, " ");

	    if ( ! value || strlen(value) == 0) {
		    cout << FG_FACE2 << "verbosity level: You must provide a value" << FG_NORMAL << endl;
		    continue;
	    }
	    else
		    verbosityLevel = atoi(value);
	    
	    if ( ! quiet )
	     cout << "New verbosity level: " << FG_FACE3 << verbosityLevel << FG_NORMAL << endl;

    }
    else if ( strstr (command, "optimize") ) {

      if (opt_optimize)
	opt_optimize = 0;
      else
	opt_optimize = 1;

      if ( ! quiet )
       cout << "Optimization: " << FG_FACE3 << ((opt_optimize) ? "yes" : "no") << FG_NORMAL << endl;

    }
    else if ( strstr (command, "affine-quad") ) {

	    if (affineQuadrics)
		    affineQuadrics = false;
	    else
		    affineQuadrics = true;

	    if ( ! quiet )
        	    cout << "Affine quadrics (in output): " << FG_FACE3 << ((affineQuadrics) ? "yes" : "no") << FG_NORMAL << endl;

    }
    else if ( strstr (command, "euclidean-type") ) {

	    if (showEuclideanType)
		    showEuclideanType = false;
	    else
		    showEuclideanType = true;

	    if ( ! quiet )
        	    cout << "Show Euclidean type of input quadrics (in output): " << FG_FACE3 << ((showEuclideanType) ? "yes" : "no") << FG_NORMAL << endl;

    }
    else if ( strstr (command, "affine-param") ) {

	    if (affineParametrizations)
		    affineParametrizations = false;
	    else
		    affineParametrizations = true;

	    if ( ! quiet )
        	    cout << "Affine parametrizations: " << FG_FACE3 << ((affineParametrizations) ? "yes" : "no") << FG_NORMAL << endl;

    }
    else if ( strstr (command, "omit") ) {

	    if (omitComponentsNotInRealSpace)
		    omitComponentsNotInRealSpace = false;
	    else
		    omitComponentsNotInRealSpace = true;

	    if ( ! quiet )
        	    cout << "Omit components which are not in R^3: " << FG_FACE3 << ((omitComponentsNotInRealSpace) ? "yes" : "no") << FG_NORMAL << endl;

    }
    else if ( strstr (command, "latex") ) {

	    if ( useLaTeX ) 
		    useLaTeX = false;
	    else 
		    useLaTeX = true;

	    if ( ! quiet )	    
		    cout << "Using LaTeX output for sqrt(delta)" << endl;
	    
    }
    else if ( strlen(command) == 0) continue;
    else {
      short i;
      if ( (i = index_of(command)) < 0 )
	cout << FG_FACE2 << "Undefined command/variable: " << command << endl << FG_NORMAL << "Check " << FG_FACE1 << "\"help\"" << FG_NORMAL << " for more information" << endl << endl;
      else {
	cout << FG_FACE3;
	print_quadric(shell_vars[i].value, cout);
	cout << FG_NORMAL;
	cout << endl;
      }
    }
  }

   #ifndef NDEBUG
     std::cout << "quitting shell_main" << std::endl;
   #endif 


}
