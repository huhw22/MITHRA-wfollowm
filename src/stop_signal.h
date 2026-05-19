#ifndef STOP_SIGNAL_H
#define STOP_SIGNAL_H

#include <csignal>

namespace MITHRA
{
  void installStopSignalHandlers();
  bool stopRequested();
}

#endif