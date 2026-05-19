#include "stop_signal.h"

namespace MITHRA
{
  namespace
  {
    volatile std::sig_atomic_t g_stopRequested = 0;

    void handleStopSignal(int)
    {
      g_stopRequested = 1;
    }
  }

  void installStopSignalHandlers()
  {
    std::signal(SIGTERM, handleStopSignal);
    std::signal(SIGINT,  handleStopSignal);
  }

  bool stopRequested()
  {
    return g_stopRequested != 0;
  }
}