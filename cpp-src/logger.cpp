#ifndef logger_cpp
#define logger_cpp
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

enum LogType { ERROR = -3, WARN = -2, OUTPUT = -1, INFO = 0, DEBUG = 1, TRACE = 2 };

class Logger {
public:
  LogType loglevel;
  ostream& of;
  string msg;
  LogType msglevel;
  
  Logger() : of (cout) {
    //    of.open("log.txt");
    //    of = cout;
    loglevel = INFO;
    msglevel = INFO;
  }
  
   Logger( LogType inlevel, ostream& os ): of( os ) {
    //of.open("log.txt");
    //    of = cout;
    loglevel = inlevel;
    msglevel = INFO;
  }

  void set_level( LogType Lin ) {
    loglevel = Lin;
  }

  void operator()( LogType level, string msg ) {
    if (level <= loglevel) {
       time_t cTime = time( NULL );
       string stamp( ctime( &cTime ) );
       stamp.pop_back();
       
      switch (level) {
      case ERROR:
	of << stamp << "\033[31m [ERROR] \033[0m" << msg << endl;
	break;
      case WARN:
	of << stamp << "\033[33m [WARN] \033[0m" << msg << endl;
	break;
      case OUTPUT:
	of << stamp << " [OUTPUT] " << msg << endl;
	break;
      case INFO:
	of << stamp << "\033[32m [INFO] \033[0m" << msg << endl;
	break;
      case DEBUG:
	of << stamp << "\033[31m [DEBUG] \033[0m" << msg << endl;
	break;
      case TRACE:
	of << stamp << "\033[31m [TRACE] \033[0m" << msg << endl;
	break;
      }
    }	     
  }

  ~Logger() {
    //of.close();
  }
};

class endlclass {
  
} endL;

template <typename T>
Logger& operator<<( Logger& lhs, T word ) {
  lhs.msg += to_string(word);
  return lhs;
}

Logger& operator<<( Logger& lhs, const char* word ) {
  string tmp( word );
  lhs.msg += tmp;
  return lhs;
}

Logger& operator<<( Logger& lhs, string word ) {
  lhs.msg += word;
  return lhs;
}

Logger& operator<<( Logger& lhs, endlclass word ) {
  lhs( lhs.msglevel, lhs.msg );
  lhs.msg.clear();

  return lhs;
}

Logger& operator<<( Logger& lhs, LogType inmsgLevel ) {
  lhs.msglevel = inmsgLevel;

  return lhs;
}



#endif
