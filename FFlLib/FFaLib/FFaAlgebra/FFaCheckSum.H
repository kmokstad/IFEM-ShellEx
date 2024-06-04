#ifndef FFA_CHECK_SUM_H
#define FFA_CHECK_SUM_H


//! \brief Dummy class doing nothing.
class FFaCheckSum
{
public:
  void reset() {}
  unsigned int getCurrent() const { return 0; }
  template<class T> void add(const T&) {}
  template<class T> void add(const T&, int) {}
};

#endif
