// This program is part of RAZIN
// Copyright Денис Крыськов 2014
// Licence: GNU General Public License (GPL)

#ifndef PROFILE__H
#define PROFILE__H

#define SHOW_TIME 1

#if SHOW_TIME
 #include <time.h>
 #define MARK_TIME(mt_x)      clock_t mt_x=clock();
 #define DUMP_TIME(mt_x,mt_y) flint_printf("%s %w\n",mt_x,clock()-mt_y);
#else
 #define MARK_TIME(x)    /**/
 #define DUMP_TIME(x,y)  /***/
#endif

#endif
