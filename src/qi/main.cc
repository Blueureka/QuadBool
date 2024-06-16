/** An interactive shell to perform live intersection calculation
 *  with the libqi library. */ 

#include <qi/parse_args.h>
#include <qi/shell.h>

int main(int argc, char	**argv)
{

  parseArgs (argc, argv);
  shell_main ();
}
