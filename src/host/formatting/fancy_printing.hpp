/**
 * @file fancy_printing.hpp
 * @author C. Noga (github.com/nohgz)
 * @brief Contains various utilities for printing fancy messages
 * @version 1.0
 * @date 2026-10-05
 */
#pragma once
#include <stdio.h>
#include <stdarg.h>

namespace host {
  namespace fancy_printing {
    #define COLOR_RED     "\033[1;31m"
    #define COLOR_GREEN   "\033[1;32m"
    #define COLOR_YELLOW  "\033[1;33m"
    #define COLOR_MAGENTA "\033[1;35m"
    #define COLOR_CYAN    "\033[1;36m"
    #define COLOR_RESET   "\033[0m"

    static inline void return_ok() {
      fprintf(stderr, "\r[ %s OK %s ]", COLOR_GREEN, COLOR_RESET);
    }

    static inline void return_fail() {
      fprintf(stderr, "\r[ %sFAIL%s ]\n", COLOR_RED, COLOR_RESET);
    }

    static inline void print_blankmsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\n[      ] ");
      vfprintf(stderr, msg, args);
      va_end(args);
    }

    static inline void okmsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\r[ %s OK %s ] ", COLOR_GREEN, COLOR_RESET);
      vfprintf(stderr, msg, args);
      fprintf(stderr, "\n");
      va_end(args);
    }

    static inline void failmsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\r[ %sFAIL%s ] ", COLOR_RED, COLOR_RESET);
      vfprintf(stderr, msg, args);
      fprintf(stderr, "\n");
      va_end(args);
    }

    static inline void timemsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\n\r[ %sTIME%s ] ", COLOR_YELLOW, COLOR_RESET);
      vfprintf(stderr, msg, args);
      va_end(args);
    }

    static inline void infomsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\n[ %sINFO%s ] ", COLOR_CYAN, COLOR_RESET);
      vfprintf(stderr, msg, args);
      va_end(args);
    }

    static inline void debugmsg(const char* msg, ...) {
      va_list args;
      va_start(args, msg);
      fprintf(stderr, "\n[ %sDBUG%s ] ", COLOR_MAGENTA, COLOR_RESET);
      vfprintf(stderr, msg, args);
      va_end(args);
    }

    static inline void header() {
      fprintf(stderr, "\n\
██████╗ ██╗  ██╗\033[38;5;160m ██████╗ \033[0m███████╗███╗   ██╗██╗██╗  ██╗\n\
██╔══██╗██║  ██║\033[38;5;196m██╔═══██╗\033[0m██╔════╝████╗  ██║██║╚██╗██╔╝\n\
██████╔╝███████║\033[38;5;202m██║   ██║\033[0m█████╗  ██╔██╗ ██║██║ ╚███╔╝ \n\
██╔═══╝ ██╔══██║\033[38;5;208m██║   ██║\033[0m██╔══╝  ██║╚██╗██║██║ ██╔██╗ \n\
██║     ██║  ██║\033[38;5;214m╚██████╔╝\033[0m███████╗██║ ╚████║██║██╔╝ ██╗\n\
╚═╝     ╚═╝  ╚═╝\033[38;5;220m ╚═════╝ \033[0m╚══════╝╚═╝  ╚═══╝╚═╝╚═╝  ╚═╝\n\n"
      );
    }
  }
}
