#ifndef _QI_ANSICOLORS_H_
#define _QI_ANSICOLORS_H_

/** Pretty printing constants */
#ifdef LIBQI_COLOR_SHELL
#define FG_BOLD    "\033[1m"
#define FG_NORMAL  "\033[0m"
#define FG_FACE1   "\033[35;1m" // Currently set to magenta + bold
#define FG_FACE2   "\033[31;1m" // Currently set to red + bold
#define FG_FACE3   "\033[34;1m" // Currently set to blue + bold
#else
#define FG_BOLD    ""
#define FG_NORMAL  ""
#define FG_FACE1   ""
#define FG_FACE2   ""
#define FG_FACE3   ""
#endif

#endif
