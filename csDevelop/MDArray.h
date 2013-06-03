/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef MDArray_h
#define MDArray_h

//
// This is a template class that provides a multi-dimensional array
// like interface that attempts to be very efficient.  This is
// essentially a proxy (wrapper) class that allows T*, T**, T***,
// etc. data be referenced from a single data type.  Reference counting
// is implemented to make sharing efficient and to prevent memory leaks.
// By choice, copy-on-write is not implemented.
//
// See Meyers, Scott "More Effective C++", Item 29 _Reference Counting_.
//
// In debug mode, this does full bounds checking as well.
//
//
#include <vector>
#include <iostream>
#include <stdexcept>


//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a vcsError is thrown. A diagnostic
 * string indicating where the error
 * occured is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param proc  Character string or std:string expression indicating the procedure
 *              where the assertion failed
 * @ingroup errorhandling
 */
#ifdef DEBUG_MODE_JJ
#define ThrowAssert(expr)  ((expr) ? (void) 0 : throw std::logic_error(std::string("MDArray: failed assert: ") + #expr))
#else
#define ThrowAssert(expr) ((void) 0);
#endif


//! Multi-dimensional Array class that allows T*, T**, and T*** data to be
//! referenced from a single data type.
/*!
 *
 *   The dimension of the array is carried within the array itself. It may be passed around without
 *   knowledge of the actual dimension. Whenever it is used, the dimension of the array must be
 *   known. 
 *
 *   In debug mode, this does full bounds checking as well.
 *
 *   Reference counting is implemented to make sharing efficient and to prevent memory leaks.
 *   By choice, copy-on-write is not implemented.
 *
 *   Internally, the data is layed out in "C" format. For example, A(i,j,k) array
 *   has the indecices k varying must rapidly within the data itself.
 *
 *   Layout of the data:
 *
 *   A(0,0,0)
 *   A(0,0,1)
 *   A(0,1,0)
 *   A(0,1,1)
 *   A(1,0,0)
 *   A(1,0,1)
 *   A(1,1,0)
 *   A(1,1,1)
 *
 */ 
template<class T>
class MDArray
{

  template<class T2>
  friend
  std::ostream & operator<<(std::ostream & os,
			    const MDArray<T2> & a);
public:
  MDArray();
  MDArray(const MDArray<T> & rhs);
  MDArray(const std::vector<int> & exts);
  virtual ~MDArray();

  //! Assignment operator for the MDArray class
  MDArray & operator=(const MDArray<T> & rhs);

  template<class T2>
  void deep_copy(const MDArray<T2> & src);
      
  static void vector_resize(std::vector< MDArray<T> > & vec,
			    const int new_size);

public:
  // 0D access.
  T & operator() ();
  const T & operator() () const;

  // 1D access.
  T & operator() (const int & i);
  const T & operator() (const int & i) const;

  // 2D access.
  T & operator() (const int & i,
		  const int & j);
  const T & operator() (const int & i,
			const int & j) const;

  // 3D access.
  T & operator() (const int & i,
		  const int & j,
		  const int & k);
  const T & operator() (const int & i,
			const int & j,
			const int & k) const;

  // 4D access.
  T & operator() (const int & i,
		  const int & j,
		  const int & k,
		  const int & l);
  const T & operator() (const int & i,
			const int & j,
			const int & k,
			const int & l) const;
  // 5D access.
  T & operator() (const int & i,
		  const int & j,
		  const int & k,
		  const int & l,
		  const int & m);
  const T & operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m) const;

  // 6D access.
  T & operator() (const int & i,
		  const int & j,
		  const int & k,
		  const int & l,
		  const int & m,
		  const int & n);
  const T & operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n) const;

  // 7D access.
  T & operator() (const int & i,
		  const int & j,
		  const int & k,
		  const int & l,
		  const int & m,
		  const int & n,
		  const int & o);
  const T & operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n,
			const int & o) const;
  // access via vector of indices
  T & operator() (const std::vector<int> & indices);
  const T & operator() (const std::vector<int> & indices) const;

public:

  //! Public const iterator class belonging to MDArray<T>
  /*!
   *
   */
  class constMDArrayIterator
  {
  public:
    constMDArrayIterator(); // disallowed

    constMDArrayIterator( const MDArray<T> & a, std::vector<int> & i )
      : indices( i ), array( a ), index0(0) {}

    constMDArrayIterator( const constMDArrayIterator & iter )
      : indices( iter.indices ), array( iter.array ), index0(iter.index0) {}

    bool operator==(const constMDArrayIterator & iter) const
    { return &array == &(iter.array) && indices == iter.indices && index0 == iter.index0; }

    bool operator!=(const constMDArrayIterator & iter) const
    { return &array != &(iter.array) || indices != iter.indices || index0 != iter.index0; }

    constMDArrayIterator & operator = ( const constMDArrayIterator & iter )
    { array = iter.array; indices = iter.indices; index0 = iter.index0; return *this ; }

    const T & operator *  () const
    { return array(indices); }

    const T * operator -> () const
    { return &(array(indices)); }

    constMDArrayIterator& operator++()
    { increment_indices(); return *this ; }

    constMDArrayIterator operator++(int)
    { constMDArrayIterator tmp(*this); ++*this ; return tmp; }

    constMDArrayIterator& operator--()
    { decrement_indices(); return *this ; }

    constMDArrayIterator operator--(int)
    { constMDArrayIterator tmp(*this); --*this ; return tmp; }

    const std::vector<int> & get_indices() const
    { return indices; }

  private:
    std::vector<int> indices;
    const MDArray<T> & array;
    int index0;

    void increment_indices()
    {
      const std::vector<int> & exts = array.get_extents();
      const int dims = array.get_dims();
      if(dims == 0)
	{
	  ++index0;
	  return;
	}
      int i = dims - 1;
      (indices[i])++;
      while ( i > 0 && indices[i] >= exts[i] )
	{
	  indices[i] = 0;
	  (indices[--i])++;
	}
    }

    void decrement_indices()
    {
      const std::vector<int> & exts = array.get_extents();
      const int dims = array.get_dims();
      if(dims == 0)
	{
	  --index0;
	  return;
	}
      int i = dims - 1;
      (indices[i])--;
      while ( i > 0 && indices[i] < 0 )
	{
	  indices[i] = exts[i] - 1;
	  (indices[--i])--;
	}
    }
  };

  class MDArrayIterator
  {
  public:
    MDArrayIterator(); // disallowed

    MDArrayIterator( MDArray<T> & a, std::vector<int> & i )
      : indices( i ), array( a ), index0( 0 ) {}

    MDArrayIterator( const MDArrayIterator & iter )
      : indices( iter.indices ), array( iter.array ), index0( iter.index0 ) {}

    bool operator==(const MDArrayIterator & iter) const
    { return &array == &(iter.array) && indices == iter.indices && index0 == iter.index0; }

    bool operator!=(const MDArrayIterator & iter) const
    { return &array != &(iter.array) || indices != iter.indices || index0 != iter.index0; }

    MDArrayIterator & operator = ( const MDArrayIterator & iter )
    { array = iter.array; indices = iter.indices; index0 = iter.index0; return *this ; }

    T & operator *  () const
    { return array(indices); }

    T * operator -> () const
    { return &(array(indices)); }

    MDArrayIterator& operator++()
    { increment_indices(); return *this ; }

    MDArrayIterator operator++(int)
    { MDArrayIterator tmp(*this); ++*this ; return tmp; }

    MDArrayIterator& operator--()
    { decrement_indices(); return *this ; }

    MDArrayIterator operator--(int)
    { MDArrayIterator tmp(*this); --*this ; return tmp; }

    const std::vector<int> & get_indices() const
    { return indices; }

  private:
    std::vector<int> indices;
    MDArray<T> & array;
    int index0;

    void increment_indices()
    {
      const std::vector<int> & exts = array.get_extents();
      const int dims = array.get_dims();
      if(dims == 0)
	{
	  ++index0;
	  return;
	}
      int i = dims - 1;
      (indices[i])++;
      while ( i > 0 && indices[i] >= exts[i] )
	{
	  indices[i] = 0;
	  (indices[--i])++;
	}
    }

    void decrement_indices()
    {
      const std::vector<int> & exts = array.get_extents();
      const int dims = array.get_dims();
      if(dims == 0)
	{
	  --index0;
	  return;
	}
      int i = dims - 1;
      (indices[i])--;
      while ( i > 0 && indices[i] < 0 )
	{
	  indices[i] = exts[i] - 1;
	  (indices[--i])--;
	}
    }
  };

  // bidirectional iterators
  typedef MDArrayIterator iterator;
  typedef constMDArrayIterator const_iterator;

  inline iterator begin();
  inline iterator end();

  inline const_iterator begin() const;
  inline const_iterator end() const;

  int get_dims() const { return data->extents.size(); }
  const std::vector<int> & get_extents() const { return data->extents; }
  const std::vector<int> & get_capacity() const { return data->capacity; }
  bool empty() { return data->empty(); }

  // These methods resize the MDArray, possibly changing its
  // dimension.  Deallocation/allocation are only done when the new
  // extents are different than the old extents so it's cheap to
  // call this method when it may not be necessary.
  void resize(const std::vector<int> & exts) { data->resize(exts); }

  template<class T2>
  void set_all(const T2 & value);
      
  //! Zero all of the entries in the array
  /*!
   * Note, zeroing is not done by default. You have to explictly zero
   * the array, or you may get undefs.
   */
  void zero_all();

private:

  //! This data class is a private class of MDArray
  /*!
   *  This class actually holds the data for the object
   *  There may be only one of these.
   */
  class MDArray_Data {
  public:
    MDArray_Data();
    MDArray_Data(const std::vector<int> & exts);
    ~MDArray_Data();

    // This actually frees the allocated memory which is NOT
    // done in the destructor.  (This permits cheap, shallow copies).
    void resize(const std::vector<int> & exts);
        
    bool empty() { return 0 == real_data; }

  private:
    void free();
    void reset_pointers();
    void allocate_pointers();
  public:
    std::vector<int> extents;
    std::vector<int> capacity;

    //! Reference counter for the data. 
    size_t ref_count;
    T * real_data;
    T * data1d;
    T ** data2d;
    T *** data3d;
    T **** data4d;
    T ***** data5d;
    T ****** data6d;
    T ******* data7d;
  };

  //! Pointer within MDArray that points to the location of the data
  //! for the array
  MDArray_Data * data;
};
// End of the MDArray description

template<class T>
void print_info(const MDArray<T> & a);

template<class T>
inline T &
MDArray<T>::operator() ()
{
  ThrowAssert(0 == data->extents.size());
  ThrowAssert(0 != data->data1d);
  return *(data->data1d);
}

template<class T>
inline const T &
MDArray<T>::operator() () const
{
  ThrowAssert(0 == data->extents.size());
  ThrowAssert(0 != data->data1d);
  return *(data->data1d);
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i)
{
  ThrowAssert(1 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(0 != data->data1d);
  return data->data1d[i];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i) const
{
  ThrowAssert(1 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(0 != data->data1d);
  return data->data1d[i];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j)
{
  ThrowAssert(2 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(0 != data->data2d);
  return data->data2d[i][j];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j) const
{
  ThrowAssert(2 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(0 != data->data2d);
  return data->data2d[i][j];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k)
{
  ThrowAssert(3 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(0 != data->data3d);
  return data->data3d[i][j][k];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k) const
{
  ThrowAssert(3 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(0 != data->data3d);
  return data->data3d[i][j][k];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l)
{
  ThrowAssert(4 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(0 != data->data4d);
  return data->data4d[i][j][k][l];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l) const
{
  ThrowAssert(4 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(0 != data->data4d);
  return data->data4d[i][j][k][l];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m)
{
  ThrowAssert(5 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(0 != data->data5d);
  return data->data5d[i][j][k][l][m];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m) const
{
  ThrowAssert(5 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(0 != data->data5d);
  return data->data5d[i][j][k][l][m];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n)
{
  ThrowAssert(6 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(n >= 0 && n < data->extents[5]);
  ThrowAssert(0 != data->data6d);
  return data->data6d[i][j][k][l][m][n];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n) const
{
  ThrowAssert(6 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(n >= 0 && n < data->extents[5]);
  ThrowAssert(0 != data->data6d);
  return data->data6d[i][j][k][l][m][n];
}

template<class T>
inline T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n,
			const int & o)
{
  ThrowAssert(7 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(n >= 0 && n < data->extents[5]);
  ThrowAssert(o >= 0 && o < data->extents[6]);
  ThrowAssert(0 != data->data7d);
  return data->data7d[i][j][k][l][m][n][o];
}

template<class T>
inline const T &
MDArray<T>::operator() (const int & i,
			const int & j,
			const int & k,
			const int & l,
			const int & m,
			const int & n,
			const int & o) const
{
  ThrowAssert(7 == data->extents.size());
  ThrowAssert(i >= 0 && i < data->extents[0]);
  ThrowAssert(j >= 0 && j < data->extents[1]);
  ThrowAssert(k >= 0 && k < data->extents[2]);
  ThrowAssert(l >= 0 && l < data->extents[3]);
  ThrowAssert(m >= 0 && m < data->extents[4]);
  ThrowAssert(n >= 0 && n < data->extents[5]);
  ThrowAssert(o >= 0 && o < data->extents[6]);
  ThrowAssert(0 != data->data7d);
  return data->data7d[i][j][k][l][m][n][o];
}

template<class T>
inline T &
MDArray<T>::operator() (const std::vector<int> & i)
{
  ThrowAssert(i.size() == (unsigned int)get_dims());

  MDArray<T> & me = *this;
  switch(get_dims())
    {
    case 0:
      return ( me() );
    case 1:
      return ( me( i[0] ) );
    case 2:
      return ( me( i[0], i[1] ) );
    case 3:
      return ( me( i[0], i[1], i[2] ) );
    case 4:
      return ( me( i[0], i[1], i[2], i[3] ) );
    case 5:
      return ( me( i[0], i[1], i[2], i[3], i[4] ) );
    case 6:
      return ( me( i[0], i[1], i[2], i[3], i[4], i[5] ) );
    case 7:
      return ( me( i[0], i[1], i[2], i[3], i[4], i[5], i[6] ) );
    default:
      throw std::logic_error("MDArray: Unsupported dimension size in MDArray: " + get_dims());
    }
}

template<class T>
inline const T &
MDArray<T>::operator() (const std::vector<int> & i) const
{
  ThrowAssert(i.size() == (unsigned int)get_dims());

  const MDArray<T> & me = *this;
  switch(get_dims())
    {
    case 0:
      return ( me() );
    case 1:
      return ( me( i[0] ) );
    case 2:
      return ( me( i[0], i[1] ) );
    case 3:
      return ( me( i[0], i[1], i[2] ) );
    case 4:
      return ( me( i[0], i[1], i[2], i[3] ) );
    case 5:
      return ( me( i[0], i[1], i[2], i[3], i[4] ) );
    case 6:
      return ( me( i[0], i[1], i[2], i[3], i[4], i[5] ) );
    case 7:
      return ( me( i[0], i[1], i[2], i[3], i[4], i[5], i[6] ) );
    default:
      throw std::logic_error("MDArray: Unsupported dimension size in MDArray: " + get_dims());
    }
}

template<class T>
inline
typename MDArray<T>::iterator
MDArray<T>::begin()
{
  std::vector<int> start( get_dims() ); // initializes to all zeros
  return MDArrayIterator( *this, start );
}

template<class T>
inline
typename MDArray<T>::iterator
MDArray<T>::end()
{
  const int dims(get_dims());
  std::vector<int> end( dims ); // initializes to all zeros
  if(dims == 0)
    {
      return ++MDArrayIterator( *this, end );
    }
  else
    {
      const std::vector<int> & exts = get_extents();
      end[0] = exts[0];
      return MDArrayIterator( *this, end );
    }
}

template<class T>
inline
typename MDArray<T>::const_iterator
MDArray<T>::begin() const
{
  std::vector<int> start( get_dims() ); // initializes to all zeros
  return constMDArrayIterator( *this, start );
}

template<class T>
inline
typename MDArray<T>::const_iterator
MDArray<T>::end() const
{
  const int dims(get_dims());
  std::vector<int> end( dims ); // initializes to all zeros
  if(dims == 0)
    {
      return ++constMDArrayIterator( *this, end );
    }
  else
    {
      const std::vector<int> & exts = get_extents();
      end[0] = exts[0];
      return constMDArrayIterator( *this, end );
    }
}

template<class T>
inline MDArray<T>::MDArray_Data::MDArray_Data() :
  ref_count(1),
  real_data(0),
  data1d(0),
  data2d(0),
  data3d(0),
  data4d(0),
  data5d(0),
  data6d(0),
  data7d(0)
{
  /* %TRACE% */  /* %TRACE% */
}

template<class T>
inline MDArray<T>::MDArray_Data::MDArray_Data(const std::vector<int> & exts) :
  ref_count(1),
  real_data(0),
  data1d(0),
  data2d(0),
  data3d(0),
  data4d(0),
  data5d(0),
  data6d(0),
  data7d(0)
{
  /* %TRACE% */  /* %TRACE% */
  resize(exts);
}

template<class T>
inline MDArray<T>::MDArray_Data::~MDArray_Data()
{
  /* %TRACE% */  /* %TRACE% */
  free();
}

template<class T>
inline void
MDArray<T>::MDArray_Data::resize(const std::vector<int> & exts)
{
  /* %TRACE% */  /* %TRACE% */

  // Allocates memory but does not reset pointer values.

  // We only need to re-allocate memory if there are a different number
  // of extents or if the new extents are larger than the current
  // capacity.  If the exts == extents then we are identical and
  // return even earlier.
  bool need_to_alloc       = false;
  bool identical           = true;
  const size_t capacity_size = capacity.size();
  const size_t ext_size      = exts.size();
  if (capacity_size != ext_size || real_data == 0) {
    identical = false;
    need_to_alloc = true;
  } else {
    for (size_t i=0; i < capacity_size; ++i) {
      if(extents[i] != exts[i]) {
	identical = false;
	if (capacity[i] < exts[i]) {
	  need_to_alloc = true;
	  break;
	}
      }
    }
  }
  
  if (identical) {
      return;
  }
        
  extents = exts;
      
  if ( !need_to_alloc) {
    return;
  }

  // Free the old data.
  free();

  // Save the new extents and capacity (free clobbers these).
  extents  = exts;
  capacity = exts;

  // Get the total amout of memory needed.
  if(capacity.size() > 7) {
    throw std::logic_error("MDArray: Error: The MDArray class is currently limited to 7 or fewer dimensions.");
  }
  
  int total_size = 1;
  for (size_t i = 0; i < capacity.size(); ++i) {
    total_size *= capacity[i];
  }
  
  // Allocate the real data as a contiguous chunk of memeory.
  try {
    real_data = new T[total_size];
  } catch (std::bad_alloc &) {
    throw std::logic_error("MDArray: MDArray: Memory allocation error");
  }

  allocate_pointers();
  reset_pointers();
}

template<class T> inline void MDArray<T>::MDArray_Data::allocate_pointers() {
  /* %TRACE% */  /* %TRACE% */
  // Allocate the multi-dimensional pointer indexes and point them into the real data array.
  try  {
    switch(capacity.size())
      {
      case 0:
      case 1:
	data1d = real_data;
	break;
      case 2:
	data2d = new T *[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data2d[i] = &real_data[capacity[1]*i];
	  }
	break;
      case 3:
	data3d = new T **[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data3d[i] = new T *[capacity[1]];
	  }
	break;
      case 4:
	data4d = new T ***[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data4d[i] = new T **[capacity[1]];
	    for(int j=0; j < capacity[1]; ++j)
	      {
		data4d[i][j] = new T *[capacity[2]];
	      }
	  }
	break;
      case 5:
	data5d = new T ****[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data5d[i] = new T ***[capacity[1]];
	    for(int j=0; j < capacity[1]; ++j)
	      {
		data5d[i][j] = new T **[capacity[2]];
		for(int k=0; k < capacity[2]; ++k)
		  {
		    data5d[i][j][k] = new T *[capacity[3]];
		  }
	      }
	  }
	break;
      case 6:
	data6d = new T *****[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data6d[i] = new T ****[capacity[1]];
	    for(int j=0; j < capacity[1]; ++j)
	      {
		data6d[i][j] = new T ***[capacity[2]];
		for(int k=0; k < capacity[2]; ++k)
		  {
		    data6d[i][j][k] = new T **[capacity[3]];
		    for(int l=0; l < capacity[3]; ++l)
		      {
			data6d[i][j][k][l] = new T *[capacity[4]];
		      }
		  }
	      }
	  }
	break;
      case 7:
	data7d = new T ******[capacity[0]];
	for(int i=0; i < capacity[0]; ++i)
	  {
	    data7d[i] = new T *****[capacity[1]];
	    for(int j=0; j < capacity[1]; ++j)
	      {
		data7d[i][j] = new T ****[capacity[2]];
		for(int k=0; k < capacity[2]; ++k)
		  {
		    data7d[i][j][k] = new T ***[capacity[3]];
		    for(int l=0; l < capacity[3]; ++l)
		      {
			data7d[i][j][k][l] = new T **[capacity[4]];
			for(int m=0; m < capacity[4]; ++m)
			  {
			    data7d[i][j][k][l][m] = new T *[capacity[5]];
			  }
		      }
		  }
	      }
	  }
	break;
      }
  } catch (std::bad_alloc &) {
    throw std::logic_error("MDArray: MDArray: Memory allocation error");
  }
}

template<class T>
inline void
MDArray<T>::MDArray_Data::reset_pointers()
{ /* %TRACE% */  /* %TRACE% */
  // Resets the pointer arrays incase the real_data space changes.
  switch(capacity.size())
    {
    case 0:
    case 1:
      data1d = real_data;
      break;
    case 2:
      for(int i=0; i < capacity[0]; ++i)
	{
	  data2d[i] = &real_data[capacity[1]*i];
	}
      break;
    case 3:
      for(int i=0; i < capacity[0]; ++i)
	{
	  for(int j=0; j < capacity[1]; ++j)
	    {
	      data3d[i][j] = &real_data[capacity[2]*(capacity[1]*i + j)];
	    }
	}
      break;
    case 4:
      for(int i=0; i < capacity[0]; ++i)
	{
	  for(int j=0; j < capacity[1]; ++j)
	    {
	      for(int k=0; k < capacity[2]; ++k)
		{
		  data4d[i][j][k] = &real_data[capacity[3]*(capacity[2]*(capacity[1]*i + j) + k)];
		}
	    }
	}
      break;
    case 5:
      for(int i=0; i < capacity[0]; ++i)
	{
	  for(int j=0; j < capacity[1]; ++j)
	    {
	      for(int k=0; k < capacity[2]; ++k)
		{
		  for(int l=0; l < capacity[3]; ++l)
		    {
		      data5d[i][j][k][l] = &real_data[capacity[4]*(capacity[3]*(capacity[2]*(capacity[1]*i + j) + k) + l)];
		    }
		}
	    }
	}
      break;
    case 6:
      for(int i=0; i < capacity[0]; ++i)
	{
	  for(int j=0; j < capacity[1]; ++j)
	    {
	      for(int k=0; k < capacity[2]; ++k)
		{
		  for(int l=0; l < capacity[3]; ++l)
		    {
		      for(int m=0; m < capacity[4]; ++m)
			{
			  data6d[i][j][k][l][m] = &real_data[capacity[5]*(capacity[4]*(capacity[3]*(capacity[2]*(capacity[1]*i + j) + k) + l) + m)];
			}
		    }
		}
	    }
	}
      break;
    case 7:
      for(int i=0; i < capacity[0]; ++i)
	{
	  for(int j=0; j < capacity[1]; ++j)
	    {
	      for(int k=0; k < capacity[2]; ++k)
		{
		  for(int l=0; l < capacity[3]; ++l)
		    {
		      for(int m=0; m < capacity[4]; ++m)
			{
			  for(int n=0; n < capacity[5]; ++n)
			    {
			      data7d[i][j][k][l][m][n] = &real_data[capacity[6]*(capacity[5]*(capacity[4]*(capacity[3]*(capacity[2]*(capacity[1]*i + j) + k) + l) + m) +n)];
			    }
			}
		    }
		}
	    }
	}
      break;
    }
}

template<class T>
inline void
MDArray<T>::MDArray_Data::free()
{ /* %TRACE% */  /* %TRACE% */
  switch(capacity.size())
    {
    case 0:
    case 1:
      // data1d is just a copy of real_data which is de-allocated below.
      break;
    case 2:
      if(0 != data2d)
	{
	  delete []data2d;
	}
      break;
    case 3:
      if(0 != data3d)
	{
	  for(int i=0; i < capacity[0]; ++i)
	    {
	      if(0 != data3d[i])
		{
		  delete []data3d[i];
		}
	    }
	  delete []data3d;
	}
      break;
    case 4:
      if(0 != data4d)
	{
	  for(int i=0; i < capacity[0]; ++i)
	    {
	      if(0 != data4d[i])
		{
		  for(int j=0; j < capacity[1]; ++j)
		    {
		      if(0 != data4d[i][j])
			{
			  delete []data4d[i][j];
			}
		    }
		  delete []data4d[i];
		}
	    }
	  delete []data4d;
	}
      break;
    case 5:
      if(0 != data5d)
	{
	  for(int i=0; i < capacity[0]; ++i)
	    {
	      if(0 != data5d[i])
		{
		  for(int j=0; j < capacity[1]; ++j)
		    {
		      if(0 != data5d[i][j])
			{
			  for(int k=0; k < capacity[2]; ++k)
			    {
			      if(0 != data5d[i][j][k])
				{
				  delete []data5d[i][j][k];
				}
			    }
			  delete []data5d[i][j];
			}
		    }
		  delete []data5d[i];
		}
	    }
	  delete []data5d;
	}
      break;
    case 6:
      if(0 != data6d)
	{
	  for(int i=0; i < capacity[0]; ++i)
	    {
	      if(0 != data6d[i])
		{
		  for(int j=0; j < capacity[1]; ++j)
		    {
		      if(0 != data6d[i][j])
			{
			  for(int k=0; k < capacity[2]; ++k)
			    {
			      if(0 != data6d[i][j][k])
				{
				  for(int l=0; l < capacity[3]; ++l)
				    {
				      if(0 != data6d[i][j][k][l])
					{
					  delete []data6d[i][j][k][l];
					}
				    }
				  delete []data6d[i][j][k];
				}
			    }
			  delete []data6d[i][j];
			}
		    }
		  delete []data6d[i];
		}
	    }
	  delete []data6d;
	}
      break;
    case 7:
      if(0 != data7d)
	{
	  for(int i=0; i < capacity[0]; ++i)
	    {
	      if(0 != data7d[i])
		{
		  for(int j=0; j < capacity[1]; ++j)
		    {
		      if(0 != data7d[i][j])
			{
			  for(int k=0; k < capacity[2]; ++k)
			    {
			      if(0 != data7d[i][j][k])
				{
				  for(int l=0; l < capacity[3]; ++l)
				    {
				      if(0 != data7d[i][j][k][l])
					{
					  for(int m=0; m < capacity[3]; ++m)
					    {
					      delete []data7d[i][j][k][l][m];
					    }
					  delete []data7d[i][j][k][l];
					}
				    }
				  delete []data7d[i][j][k];
				}
			    }
			  delete []data7d[i][j];
			}
		    }
		  delete []data7d[i];
		}
	    }
	  delete []data7d;
	}
      break;
    default:
      throw std::logic_error("MDArray: I can't be here!");
    }
  if(0 != real_data)
    {
      // FIXME: valgrind complains about this.
      delete []real_data;
    }
  extents.clear();
  capacity.clear();
  real_data = 0;
  data1d = 0;
  data2d = 0;
  data3d = 0;
  data4d = 0;
  data5d = 0;
  data6d = 0;
  data7d = 0;
}

// Wrapper Methods;
// New ctor
template<class T>
MDArray<T>::MDArray():
  data(new MDArray_Data)
{
  /* %TRACE% */  /* %TRACE% */
}

template<class T>
MDArray<T>::MDArray(const std::vector<int> & exts) :
  data(new MDArray_Data(exts))
{
  /* %TRACE% */  /* %TRACE% */
}

// Copy ctor
template<class T>
MDArray<T>::MDArray(const MDArray<T> & rhs) :
  data(rhs.data)
{
  /* %TRACE% */  /* %TRACE% */
  ++data->ref_count;
}

// Assignment operator for the MDArray 
template<class T>
MDArray<T> &
MDArray<T>::operator=(const MDArray<T> & rhs)
{
  if (data == rhs.data)  {
    return *this;
  }
  // Decrement our data reference count.
  // Looking at our data, if we are the only copy,
  // then we will not need the data any more, so delete it
  if (--data->ref_count == 0) {
    delete data;
  }
  // Copy the other guys data -> just a shallow pointer copy
  data = rhs.data;
  // Increase the ref count on the other guy's data
  ++data->ref_count;
  return *this;
}

template<class T>
MDArray<T>::~MDArray()
{
  /* %TRACE% */  /* %TRACE% */
  if(--data->ref_count == 0)
    {
      delete data;
    }
}

template<class T2>
std::ostream & operator<<(std::ostream & os, const MDArray<T2> & a)
{
  std::vector<int> exts = a.get_extents();
  switch(a.get_dims())
    {
    case 0:
      os << " = " << a() << std::endl;
      break;
    case 1:
      for(int i=0; i < exts[0]; ++i)
	os << "[" << i << "] = " << a(i) << std::endl;
      break;
    case 2:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  os << "[" << i << "]"
	     << "[" << j << "]"
	     << " = " << a(i,j) << std::endl;
      break;
    case 3:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    os << "[" << i << "]"
	       << "[" << j << "]"
	       << "[" << k << "]"
	       << " = " << a(i,j,k) << std::endl;
      break;
    case 4:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      os << "[" << i << "]"
		 << "[" << j << "]"
		 << "[" << k << "]"
		 << "[" << l << "]"
		 << " = " << a(i,j,k,l) << std::endl;
      break;
    case 5:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		os << "[" << i << "]"
		   << "[" << j << "]"
		   << "[" << k << "]"
		   << "[" << l << "]"
		   << "[" << m << "]"
		   << " = " << a(i,j,k,l,m) << std::endl;
      break;
    case 6:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		for(int n=0; n < exts[5]; ++n)
		  os << "[" << i << "]"
		     << "[" << j << "]"
		     << "[" << k << "]"
		     << "[" << l << "]"
		     << "[" << m << "]"
		     << "[" << n << "]"
		     << " = " << a(i,j,k,l,m,n) << std::endl;
      break;
    case 7:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		for(int n=0; n < exts[5]; ++n)
		  for(int o=0; o < exts[6]; ++o)
		    os << "[" << i << "]"
		       << "[" << j << "]"
		       << "[" << k << "]"
		       << "[" << l << "]"
		       << "[" << m << "]"
		       << "[" << n << "]"
		       << "[" << o << "]"
		       << " = " << a(i,j,k,l,m,n,o) << std::endl;
      break;
    default:
      throw std::logic_error("MDArray: I can't be here!");
    }
  return os;
}

template<class T>
template<class T2>
void
MDArray<T>::set_all(const T2 & value)
{

  // Determine the total capcity.
  int total_size  = 1;
  const int dims  = get_dims();
  const int * cap = &(data->capacity[0]);
  for(int i=0; i < dims; ++i)
    {
      total_size *= *cap;
      ++cap;
    }

  // Set the value.
  T * dat = &(data->real_data[0]);
  const int & const_total_size = total_size;
  for(int i=0; i < const_total_size; ++i)
    {
      *dat = value;
      ++dat;
    }
}
    
template<class T>
void
MDArray<T>::zero_all()
{

  // Determine the total capcity.
  int total_size  = 1;
  const int dims  = get_dims();
  const int * cap = &(data->capacity[0]);
  for(int i=0; i < dims; ++i)
    {
      total_size *= *cap;
      ++cap;
    }

  // Set the value.
  T * dat = &(data->real_data[0]);
  std::memset(dat,0,total_size*sizeof(T));
}

template<class T>
template<class T2>
void
MDArray<T>::deep_copy(const MDArray<T2> & src)
{
  const std::vector<int> exts = src.get_extents();
  resize( exts );

  MDArray<T> & dest = *this;
  switch(get_dims())
    {
    case 0:
      dest() = src();
      break;
    case 1:
      for(int i=0; i < exts[0]; ++i)
	dest(i) = src(i);
      break;
    case 2:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  dest(i,j) = src(i,j);
      break;
    case 3:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    dest(i,j,k) = src(i,j,k);
      break;
    case 4:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      dest(i,j,k,l) = src(i,j,k,l);
      break;
    case 5:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		dest(i,j,k,l,m) = src(i,j,k,l,m);
      break;
    case 6:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		for(int n=0; n < exts[5]; ++n)
		  dest(i,j,k,l,m,n) = src(i,j,k,l,m,n);
      break;
    case 7:
      for(int i=0; i < exts[0]; ++i)
	for(int j=0; j < exts[1]; ++j)
	  for(int k=0; k < exts[2]; ++k)
	    for(int l=0; l < exts[3]; ++l)
	      for(int m=0; m < exts[4]; ++m)
		for(int n=0; n < exts[5]; ++n)
		  for(int o=0; o < exts[6]; ++o)
		    dest(i,j,k,l,m,n,o) = src(i,j,k,l,m,n,o);
      break;
    default:
      throw std::logic_error("MDArray: Unsupported dimension size in MDArray: " + get_dims());
    }
}

template<class T>
void print_info(const MDArray<T> & a)
{
  int dims = a.get_dims();
  std::cout << "The passed in arrays has " << dims << " dimensions." << std::endl;
  std::vector<int> exts = a.get_extents();
  std::vector<int> caps = a.get_capacity();
  for(int i=0; i < dims; ++i)
    {
      std::cout << "\tdimension "  << i
		<< ": extent = "   << exts[i]
		<< ", capacity = " << caps[i]
		<< std::endl;
    }
}
 
template<class T>
void
MDArray<T>::vector_resize(std::vector< MDArray<T> > & vec,
			  const int new_size)
{
  const int old_size = vec.size();
  if ( new_size < old_size )
    {
      vec.resize(new_size);
    }
  else if ( new_size > old_size )
    {
      vec.reserve(new_size);
      for (int i=old_size; i<new_size; ++i)
	{
	  MDArray<T> new_array;
	  vec.push_back(new_array);
	}
    }
}

typedef MDArray<double> Real_MDArray;
typedef MDArray<int> Int_MDArray;



#endif // MDArray_h
